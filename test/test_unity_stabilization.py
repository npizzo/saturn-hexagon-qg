#!/usr/bin/env python3
"""Offline checks of the new campaign wrapper; no Unity or network access."""

import hashlib
import json
from pathlib import Path
import shutil
import tarfile
import unittest

from unity_slurm_tests import SlurmHarness


class StabilizationTests(SlurmHarness):
    def setUp(self):
        super().setUp()
        for key in ('SLURM_JOB_ID', 'SLURM_JOB_ACCOUNT', 'SLURM_CPUS_PER_TASK'):
            self.env.pop(key)
        previous = self.root / 'saturn-hexagon-qg'
        previous.mkdir()
        shutil.move(str(self.unity_root), str(previous / '.unity'))
        self.unity_root = previous / '.unity'
        self.unity_julia = self.unity_root / 'runtime/julia-1.10.12/bin/julia'
        with self.unity_julia.open('a') as stream:
            stream.write('''
if name == "julia" and any(value.endswith("/validate_stabilization.jl") for value in args):
    output = Path(args[args.index("--output") + 1])
    output.mkdir()
    (output / "COMPLETED").write_text("mock ordinary run completed\\n")
    (output / "CPU_GPU_VALIDATION_PASSED").write_text("mock comparison passed\\n")
    (output / "cpu_gpu_validation.toml").write_text("passed = true\\n")
''')
        self.suite = self.project / 'configs/stabilization_01'
        self.suite.mkdir()
        (self.project / 'src').mkdir()
        (self.project / 'src/Model.jl').write_text('# production fixture\n')
        for name in ('make_stabilization.jl', 'make_experiments.jl', 'validate_stabilization.jl'):
            (self.project / 'scripts' / name).write_text('# generator fixture\n')
        self.provenance_files = ['Project.toml', 'Manifest.toml', 'src/Model.jl',
            'scripts/make_stabilization.jl', 'scripts/make_experiments.jl', 'scripts/run.jl',
            'scripts/validate_stabilization.jl', 'slurm/unity/stabilization-validation.sbatch',
            'slurm/common.sh', 'slurm/ensemble.sbatch', 'slurm/unity/environment.sh',
            'slurm/unity/submit.sh', 'slurm/unity/ensemble.sbatch',
            'slurm/unity/experiments_01.sh', 'slurm/unity/experiments_01.py',
            'slurm/unity/stabilization.sh', 'slurm/unity/stabilization.py']
        self.write_provenance()
        (self.suite / 'campaign.toml').write_text('name = "offline stabilization fixture"\n')
        (self.suite / 'physics-audit.tsv').write_text('case\tinterpretation\nfixture\ttest only\n')
        (self.suite / 'stages.tsv').write_text('validation\t3\tGPU smoke checks before science\ncontrols\t2\tReview validation first\n')
        for stage, count in (('validation', 3), ('controls', 2)):
            directory = self.suite / stage
            directory.mkdir()
            members = []
            for index in range(count):
                config = directory / ('member-{}.toml'.format(index + 1))
                config.write_text('[experiment]\nname = "{}-{}"\n'.format(stage, index + 1))
                members.append(self.sha(config) + '\t' + config.name + '\n')
            (directory / 'members.tsv').write_text(''.join(members))

    @staticmethod
    def sha(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()

    def write_provenance(self):
        (self.suite / 'provenance.tsv').write_text(''.join(self.sha(self.project / name) + '\t' + name + '\n' for name in self.provenance_files))

    def call(self, command, *options, stage='validation', **kwargs):
        return self.run_script('slurm/unity/stabilization.sh', command, self.suite, stage, *options, **kwargs)

    def attempt(self, stage='validation'):
        pointers = list((self.unity_root / 'experiments-01').glob('*/' + stage + '/latest'))
        self.assertEqual(len(pointers), 1)
        return pointers[0].parent / pointers[0].read_text().strip()

    def task_environment(self):
        environment = self.calls('sbatch')[-1]['env'].copy()
        environment.update({'SLURM_JOB_ID': '12345', 'SLURM_JOB_ACCOUNT': 'unity-default',
                            'SLURM_CPUS_PER_TASK': '8', 'SLURM_ARRAY_JOB_ID': '12345', 'SLURM_ARRAY_TASK_ID': '1'})
        return environment

    def test_validate_checks_without_preparing_or_submitting(self):
        result = self.call('validate')
        self.assertIn('3 configuration hashes', result.stdout)
        self.assertFalse((self.unity_root / 'experiments-01').exists())
        self.assertFalse(self.calls('sbatch'))
        self.assertFalse(self.calls('julia'))

    def test_submission_reuses_old_runtime_but_separates_outputs(self):
        self.call('submit', env={'SATURN_UNITY_ROOT': '/stale/delta/runtime',
            'SATURN_UNITY_RUN_ROOT': '/stale/campaign/results', 'SATURN_RUN_ROOT': '/other/results',
            'JULIA_BIN': '/wrong/julia', 'SATURN_UNITY_ACCOUNT': 'explicit-account', 'SBATCH_ACCOUNT': 'wrong-account'})
        call = self.calls('sbatch')[0]
        self.assertEqual(len(self.calls('sbatch')), 1)
        self.assertIn('--array=1-3%2', call['args'])
        self.assertIn('--account=explicit-account', call['args'])
        self.assertIn('--partition=gpu-highvram,uri-gpu', call['args'])
        self.assertEqual(call['args'][-1], str(self.project / 'slurm/unity/stabilization-validation.sbatch'))
        self.assertEqual(call['env']['JULIA_BIN'], str(self.unity_julia))
        self.assertEqual(call['env']['SATURN_RUN_ROOT'], str(self.project / 'runs/stabilization-01'))
        self.assertEqual(call['env']['SATURN_PROJECT'], str(self.project))
        self.assertNotIn('SBATCH_ACCOUNT', call['env'])
        self.assertTrue((self.attempt() / 'frozen/physics-audit.tsv').is_file())
        self.assertFalse(list((self.unity_root / 'experiments-01').glob('*/controls')))

    def test_new_generator_and_wrapper_must_be_in_provenance(self):
        for name in ('scripts/make_stabilization.jl', 'slurm/unity/stabilization.sh', 'scripts/validate_stabilization.jl'):
            self.provenance_files.remove(name)
            self.write_provenance()
            self.call('submit', expected='failure')
            self.provenance_files.append(name)
        self.assertFalse(self.calls('sbatch'))

    def test_member_and_generator_tampering_fail_before_submission(self):
        for path in (self.project / 'scripts/make_stabilization.jl', self.suite / 'validation/member-1.toml'):
            original = path.read_bytes()
            path.write_bytes(original + b'# changed\n')
            self.call('submit', expected='failure')
            path.write_bytes(original)
        self.assertFalse(self.calls('sbatch'))

    def test_validation_task_uses_comparison_launcher_and_frozen_hash(self):
        self.call('submit')
        self.run_script('slurm/unity/stabilization-validation.sbatch', env=self.task_environment())
        launches = [r for r in self.calls('julia') if str(self.project / 'scripts/validate_stabilization.jl') in r['args']]
        self.assertEqual(len(launches), 1)
        self.assertFalse([r for r in self.calls('julia') if str(self.project / 'scripts/run.jl') in r['args']])
        self.assertIn(str(self.project / 'runs/stabilization-01'), ' '.join(launches[0]['args']))
        frozen = self.attempt() / 'frozen/configs/member-1.toml'
        frozen.chmod(0o644)
        frozen.write_text('# tampered\n')
        before = len(self.calls('srun'))
        self.run_script('slurm/unity/stabilization-validation.sbatch', env=self.task_environment(), expected='failure')
        self.assertEqual(len(self.calls('srun')), before)

    def test_science_stage_retains_existing_gpu_launcher(self):
        self.call('submit', stage='controls')
        self.assertEqual(self.calls('sbatch')[0]['args'][-1], str(self.project / 'slurm/unity/ensemble.sbatch'))
        self.run_script('slurm/unity/ensemble.sbatch', env=self.task_environment())
        self.assert_gpu_launches()

    def test_failed_comparison_exits_nonzero_and_collection_refuses(self):
        self.call('submit')
        environment = self.task_environment()
        environment.update({'MOCK_FAIL_TOOL': 'julia', 'MOCK_FAIL_MATCH': 'validate_stabilization.jl'})
        self.run_script('slurm/unity/stabilization-validation.sbatch', env=environment, expected=37)
        self.call('collect', expected='failure')
        self.assertFalse(list(self.attempt().glob('results-*.tar.gz')))

    def test_zero_exit_without_success_markers_is_rejected(self):
        self.call('submit')
        shutil.copy2(self.bin / 'julia', self.unity_julia)
        result = self.run_script('slurm/unity/stabilization-validation.sbatch', env=self.task_environment(), expected='failure')
        self.assertIn('both success markers', result.stderr)
        self.call('collect', expected='failure')

    def test_validation_results_collect_with_frozen_inputs(self):
        self.call('submit')
        for task in (1, 2, 3):
            environment = self.task_environment()
            environment['SLURM_ARRAY_TASK_ID'] = str(task)
            self.run_script('slurm/unity/stabilization-validation.sbatch', env=environment)
        self.call('collect')
        with tarfile.open(self.attempt() / 'results-12345.tar.gz') as archive:
            names = archive.getnames()
        self.assertIn('inputs/physics-audit.tsv', names)
        self.assertEqual(sum(name.endswith('/CPU_GPU_VALIDATION_PASSED') for name in names), 3)
        self.call('collect', expected='failure')

    def test_collection_requires_comparison_marker_not_only_completed_run(self):
        self.call('submit')
        for task in (1, 2, 3):
            environment = self.task_environment()
            environment['SLURM_ARRAY_TASK_ID'] = str(task)
            self.run_script('slurm/unity/stabilization-validation.sbatch', env=environment)
        marker = next((self.project / 'runs/stabilization-01').glob('array-*/results/CPU_GPU_VALIDATION_PASSED'))
        marker.unlink()
        self.call('collect', expected='failure')
        self.assertFalse(list(self.attempt().glob('results-*.tar.gz')))

    def test_repeat_resume_and_array_override_are_rejected(self):
        self.call('submit')
        self.call('submit', expected='failure')
        self.call('submit', '--new-attempt', env={'SATURN_RESUME': '/unrelated/checkpoint'}, expected='failure')
        self.call('submit', '--array=1-99', expected='failure')
        self.assertEqual(len(self.calls('sbatch')), 1)

    def test_failure_retains_reservation_and_reports_nonzero(self):
        self.call('submit', env={'MOCK_FAIL_TOOL': 'sbatch'}, expected='failure')
        self.assertTrue((self.attempt() / 'submission-started.json').is_file())
        self.assertFalse((self.attempt() / 'receipt.json').exists())
        self.call('submit', expected='failure')
        self.assertEqual(len(self.calls('sbatch')), 1)

    def test_status_uses_receipt_after_source_edit_and_ignores_shared_pointer(self):
        self.call('submit')
        (self.unity_root / 'last-job-id').write_text('999999\n')
        (self.project / 'scripts/make_stabilization.jl').write_text('# later edit\n')
        self.call('status')
        for tool in ('squeue', 'sacct'):
            self.assertIn('12345', self.calls(tool)[0]['args'])
            self.assertNotIn('999999', self.calls(tool)[0]['args'])

    def test_missing_runtime_and_original_checkout_use_are_rejected(self):
        self.unity_julia.unlink()
        self.call('submit', expected='failure')
        self.call('validate', env={'SATURN_STABILIZATION_RUNTIME_ROOT': self.project / '.unity'}, expected='failure')
        self.assertFalse(self.calls('sbatch'))


if __name__ == '__main__':
    unittest.main(verbosity=2)
