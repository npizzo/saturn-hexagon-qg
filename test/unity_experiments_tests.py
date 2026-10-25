#!/usr/bin/env python3
"""Offline checks of immutable experiment-stage preparation and submission.

Uses the existing isolated mock scheduler harness; never accesses Unity.
"""

import hashlib
import json
from pathlib import Path
import tarfile
import unittest

from unity_slurm_tests import SlurmHarness


class ExperimentTests(SlurmHarness):
    def setUp(self):
        super().setUp()
        for key in ('SLURM_JOB_ID', 'SLURM_JOB_ACCOUNT', 'SLURM_CPUS_PER_TASK'):
            self.env.pop(key)
        self.suite = self.project / 'configs' / 'suite with spaces'
        self.suite.mkdir()
        (self.project / 'src').mkdir()
        (self.project / 'src/Model.jl').write_text('# source fixture\n')
        (self.project / 'scripts/make_experiments.jl').write_text('# generator fixture\n')
        files = ['src/Model.jl', 'Project.toml', 'Manifest.toml', 'scripts/make_experiments.jl',
                 'slurm/unity/experiments_01.py', 'slurm/unity/ensemble.sbatch']
        (self.suite / 'provenance.tsv').write_text(''.join(self.sha(self.project / name) + '\t' + name + '\n' for name in files))
        (self.suite / 'campaign.toml').write_text('name = "offline experiment suite"\n')
        (self.suite / 'stages.tsv').write_text('screen\t2\tReview screen before choosing another stage.\nexpanded\t1\tExplicit later choice.\n')
        for stage, count in (('screen', 2), ('expanded', 1)):
            directory = self.suite / stage
            directory.mkdir()
            members = []
            for index in range(count):
                config = directory / ('member {}.toml'.format(index + 1))
                config.write_text('[experiment]\nname = "{}-{}"\n'.format(stage, index + 1))
                members.append(self.sha(config) + '\t' + config.name + '\n')
            (directory / 'members.tsv').write_text(''.join(members))

    @staticmethod
    def sha(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()

    def call(self, command, *options, stage='screen', **kwargs):
        return self.run_script('slurm/unity/experiments_01.sh', command, self.suite, stage, *options, **kwargs)

    def attempt(self, stage='screen'):
        matches = list((self.unity_root / 'experiments-01').glob('*/' + stage + '/latest'))
        self.assertEqual(len(matches), 1)
        return matches[0].parent / matches[0].read_text().strip()

    def task_environment(self):
        environment = self.calls('sbatch')[-1]['env'].copy()
        environment.update({'SLURM_JOB_ID': '12345', 'SLURM_JOB_ACCOUNT': 'unity-default',
                            'SLURM_CPUS_PER_TASK': '8', 'SLURM_ARRAY_JOB_ID': '12345', 'SLURM_ARRAY_TASK_ID': '1'})
        return environment

    def make_completed(self):
        run_root = self.project / 'runs/unity'
        for index in (1, 2):
            output = run_root / ('array-12345-{}.fixture'.format(index))
            (output / 'results').mkdir(parents=True)
            (output / 'results/COMPLETED').write_text('completed\n')
            (output / 'config.toml').write_bytes((self.attempt() / 'frozen/configs' / ('member {}.toml'.format(index))).read_bytes())
            (run_root / 'logs' / ('saturn-ensemble-12345_{}.out'.format(index))).write_text('fixture log\n')

    def test_prepare_is_local_freezes_absolute_configs_and_reuses(self):
        self.call('prepare')
        attempt = self.attempt()
        rows = (attempt / 'frozen/config-list.tsv').read_text().splitlines()
        self.assertEqual(len(rows), 2)
        self.assertTrue(all(Path(row.split('\t')[1]).is_absolute() for row in rows))
        self.assertTrue(all(str(attempt / 'frozen/configs') in row for row in rows))
        self.assertEqual((attempt / 'frozen/configs/member 1.toml').stat().st_mode & 0o222, 0)
        self.call('prepare')
        self.assertEqual(self.attempt(), attempt)
        self.assertFalse(self.calls('sbatch'))
        self.assertFalse(self.calls('julia'))

    def test_submit_one_stage_default_limit_preserves_scheduler_settings(self):
        self.call('submit', env={'SATURN_UNITY_ACCOUNT': 'verified-account', 'SBATCH_ARRAY_INX': '90-100'})
        calls = self.calls('sbatch')
        self.assertEqual(len(calls), 1)
        args = calls[0]['args']
        self.assertIn('--array=1-2%2', args)
        self.assertIn('--time=01:00:00', args)
        self.assertIn('--account=verified-account', args)
        self.assertIn('--partition=gpu-highvram,uri-gpu', args)
        self.assertIn('--constraint=x86_64&a100', args)
        self.assertFalse(any(k.startswith('SBATCH_') for k in calls[0]['env']))
        receipt = json.loads((self.attempt() / 'receipt.json').read_text())
        self.assertEqual(receipt['job_id'], '12345')
        self.assertEqual(receipt['stage'], 'screen')
        self.assertFalse(list((self.unity_root / 'experiments-01').glob('*/expanded')))

    def test_explicit_options_and_new_attempt(self):
        self.call('submit', '--concurrency', '1', '--time', '00:30:00')
        first = self.attempt()
        self.call('submit', '--new-attempt', env={'MOCK_SBATCH_SEQUENCE': '1'})
        self.assertNotEqual(self.attempt(), first)
        self.assertEqual(len(self.calls('sbatch')), 2)
        self.assertIn('--array=1-2%1', self.calls('sbatch')[0]['args'])
        self.assertEqual(self.calls('sbatch')[0]['args'].count('--time=00:30:00'), 1)
        self.assertEqual(json.loads((first / 'receipt.json').read_text())['job_id'], '12345')

    def test_repeat_and_resume_are_rejected_without_extra_jobs(self):
        self.call('submit')
        self.call('submit', expected='failure')
        self.call('submit', '--new-attempt', env={'SATURN_RESUME': '/unrelated/checkpoint'}, expected='failure')
        self.assertEqual(len(self.calls('sbatch')), 1)

    def test_submission_failure_retains_reservation_and_no_blind_retry(self):
        self.call('submit', env={'MOCK_FAIL_TOOL': 'sbatch'}, expected='failure')
        attempt = self.attempt()
        self.assertTrue((attempt / 'submission-started.json').is_file())
        self.assertFalse((attempt / 'receipt.json').exists())
        self.assertIn('intentional mock failure', (attempt / 'submission.log').read_text())
        self.call('submit', expected='failure')
        result = self.call('status')
        self.assertIn('No confirmed job receipt', result.stdout)
        self.assertEqual(len(self.calls('sbatch')), 1)
        self.assertFalse(self.calls('squeue'))

    def test_current_source_and_dependency_hash_mismatches_prevent_submit(self):
        for name in ('src/Model.jl', 'Manifest.toml'):
            path = self.project / name
            original = path.read_bytes()
            path.write_bytes(original + b'# changed\n')
            self.call('submit', expected='failure')
            path.write_bytes(original)
        self.assertFalse(self.calls('sbatch'))

    def test_changed_member_or_count_and_path_traversal_prevent_submit(self):
        manifest = self.suite / 'screen/members.tsv'
        original = manifest.read_text()
        for changed in (original.replace('member 1.toml', '../member 1.toml'), original.splitlines()[0] + '\n', original.replace(original[:64], '0' * 64, 1)):
            manifest.write_text(changed)
            self.call('submit', expected='failure')
        manifest.write_text(original)
        self.assertFalse(self.calls('sbatch'))

    def test_changed_suite_after_prepare_requires_explicit_attempt(self):
        self.call('prepare')
        (self.suite / 'campaign.toml').write_text('name = "changed"\n')
        self.call('submit', expected='failure')
        self.assertFalse(self.calls('sbatch'))

    def test_array_override_and_invalid_stages_are_not_forwarded(self):
        self.call('submit', '--array=1-99', expected='failure')
        self.call('submit', '--export=NONE', expected='failure')
        self.call('submit', stage='../screen', expected='failure')
        self.call('submit', stage='missing', expected='failure')
        self.assertFalse(self.calls('sbatch'))

    def test_per_task_verifies_current_source_before_julia(self):
        self.call('submit')
        environment = self.task_environment()
        (self.project / 'src/Model.jl').write_text('# changed after submission\n')
        self.run_script('slurm/unity/ensemble.sbatch', env=environment, expected='failure')
        self.assertFalse(self.calls('srun'))
        self.assertFalse(self.calls('julia'))

    def test_per_task_rejects_changed_frozen_config(self):
        self.call('submit')
        environment = self.task_environment()
        path = self.attempt() / 'frozen/configs/member 1.toml'
        path.chmod(0o644)
        path.write_text('# tampered\n')
        self.run_script('slurm/unity/ensemble.sbatch', env=environment, expected='failure')
        self.assertFalse(self.calls('srun'))

    def test_per_task_validates_then_uses_existing_gpu_launcher(self):
        self.call('submit')
        self.run_script('slurm/unity/ensemble.sbatch', env=self.task_environment())
        self.assert_gpu_launches()

    def test_status_uses_receipt_not_last_job_and_queries_array_elements(self):
        self.call('submit')
        (self.unity_root / 'last-job-id').write_text('999999\n')
        before = sorted(p.relative_to(self.attempt()).as_posix() for p in self.attempt().rglob('*'))
        self.call('status')
        after = sorted(p.relative_to(self.attempt()).as_posix() for p in self.attempt().rglob('*'))
        self.assertEqual(before, after)
        for tool in ('squeue', 'sacct'):
            args = self.calls(tool)[0]['args']
            self.assertIn('12345', args)
            self.assertNotIn('999999', args)
            self.assertNotIn('-X', args)
        self.assertEqual(len(self.calls('sbatch')), 1)

    def test_collect_requires_complete_members_and_no_overwrite(self):
        self.call('submit')
        self.call('collect', expected='failure')
        self.make_completed()
        self.call('collect')
        archive = self.attempt() / 'results-12345.tar.gz'
        with tarfile.open(archive) as stream:
            names = stream.getnames()
            self.assertIn('inputs/configs/member 1.toml', names)
            self.assertIn('runs/array-12345-2.fixture/results/COMPLETED', names)
            mapping = json.load(stream.extractfile('collection.json'))
            self.assertTrue(mapping['all_members_completed'])
            self.assertEqual(mapping['members'][1]['config'], 'member 2.toml')
        original = self.sha(archive)
        self.call('collect', expected='failure')
        self.assertEqual(self.sha(archive), original)
        self.assertEqual(len(self.calls('sbatch')), 1)

    def test_collection_checks_staged_config_but_allows_later_source_edits(self):
        self.call('submit')
        self.make_completed()
        staged = self.project / 'runs/unity/array-12345-1.fixture/config.toml'
        original = staged.read_bytes()
        staged.write_text('# unrelated run\n')
        self.call('collect', expected='failure')
        staged.write_bytes(original)
        (self.project / 'src/Model.jl').write_text('# later development does not invalidate saved inputs\n')
        self.call('collect')

    def test_retrieval_rejects_modified_frozen_inputs(self):
        self.call('submit')
        self.make_completed()
        config = self.attempt() / 'frozen/configs/member 1.toml'
        config.chmod(0o644)
        config.write_text('# tampered frozen input\n')
        self.call('status', expected='failure')
        self.call('collect', expected='failure')
        self.assertFalse(self.calls('squeue'))

    def test_optional_physics_audit_is_frozen_hashed_and_collected(self):
        audit = self.suite / 'physics-audit.tsv'
        audit.write_text('physics\tinterpretation\nsmall_ld\tMathematical stress test\n')
        self.call('submit')
        frozen = self.attempt() / 'frozen'
        self.assertEqual((frozen / audit.name).read_bytes(), audit.read_bytes())
        self.assertIn(self.sha(audit) + '\t' + audit.name + '\n', (frozen / 'integrity.tsv').read_text())
        self.assertEqual((frozen / audit.name).stat().st_mode & 0o222, 0)
        self.make_completed()
        self.call('collect')
        with tarfile.open(self.attempt() / 'results-12345.tar.gz') as stream:
            self.assertEqual(stream.extractfile('inputs/physics-audit.tsv').read(), audit.read_bytes())
        (frozen / audit.name).chmod(0o644)
        (frozen / audit.name).write_text('# changed interpretation\n')
        self.call('status', expected='failure')

    def test_optional_audit_addition_removal_and_changes_block_prepared_reuse(self):
        self.call('prepare')
        audit = self.suite / 'physics-audit.tsv'
        original = 'physics\tinterpretation\nsmall_ld\tMathematical stress test\n'
        audit.write_text(original)
        self.call('submit', expected='failure')  # Newly added since preparation.
        self.call('prepare', '--new-attempt')
        audit.unlink()
        self.call('submit', expected='failure')  # Removed since preparation.
        audit.write_text(original + 'another\tUpdated annotation\n')
        self.call('submit', expected='failure')  # Changed since preparation.
        audit.write_text(original)
        self.call('prepare')
        self.assertFalse(self.calls('sbatch'))


if __name__ == '__main__':
    unittest.main(verbosity=2)
