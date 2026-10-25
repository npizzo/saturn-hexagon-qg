#!/usr/bin/env python3
"""Exercise Rostami's wrapper through the existing offline Slurm harness."""
import shutil
import unittest
from test_unity_stabilization import StabilizationTests


class RostamiTests(StabilizationTests):
    def setUp(self):
        super().setUp()
        original = self.suite
        self.suite = self.project / 'configs/rostami_01'
        shutil.move(str(original), str(self.suite))
        (self.project / 'scripts/make_rostami.jl').write_text('# new generator fixture\n')
        self.provenance_files.extend(('scripts/make_rostami.jl', 'slurm/unity/rostami.sh', 'slurm/unity/rostami.py'))
        self.write_provenance()

    def call(self, command, *options, stage='validation', **kwargs):
        return self.run_script('slurm/unity/rostami.sh', command, self.suite, stage, *options, **kwargs)

    def test_submission_reuses_old_runtime_but_separates_outputs(self):
        self.call('submit', env={'SATURN_UNITY_ROOT': '/stale/delta/runtime',
            'SATURN_STABILIZATION_RUN_ROOT': '/stale/stabilization', 'SATURN_UNITY_RUN_ROOT': '/stale/results',
            'SATURN_RUN_ROOT': '/other/results', 'JULIA_BIN': '/wrong/julia', 'SBATCH_ACCOUNT': 'wrong-account'})
        call = self.calls('sbatch')[0]
        self.assertEqual(len(self.calls('sbatch')), 1)
        self.assertIn('--array=1-3%2', call['args'])
        self.assertEqual(call['args'][-1], str(self.project / 'slurm/unity/stabilization-validation.sbatch'))
        self.assertEqual(call['env']['JULIA_BIN'], str(self.unity_julia))
        self.assertEqual(call['env']['SATURN_RUN_ROOT'], str(self.project / 'runs/rostami-01'))
        self.assertNotIn('SBATCH_ACCOUNT', call['env'])
        self.assertFalse(list((self.unity_root / 'experiments-01').glob('*/controls')))

    def test_new_generator_and_wrapper_must_be_in_provenance(self):
        for name in ('scripts/make_rostami.jl', 'slurm/unity/rostami.sh', 'slurm/unity/rostami.py'):
            self.provenance_files.remove(name)
            self.write_provenance()
            self.call('submit', expected='failure')
            self.provenance_files.append(name)
        self.assertFalse(self.calls('sbatch'))

    def test_validation_task_uses_comparison_launcher_and_frozen_hash(self):
        self.call('submit')
        self.run_script('slurm/unity/stabilization-validation.sbatch', env=self.task_environment())
        launches = [r for r in self.calls('julia') if str(self.project / 'scripts/validate_stabilization.jl') in r['args']]
        self.assertEqual(len(launches), 1)
        self.assertIn(str(self.project / 'runs/rostami-01'), ' '.join(launches[0]['args']))
        frozen = self.attempt() / 'frozen/configs/member-1.toml'
        frozen.chmod(0o644)
        frozen.write_text('# changed\n')
        before = len(self.calls('srun'))
        self.run_script('slurm/unity/stabilization-validation.sbatch', env=self.task_environment(), expected='failure')
        self.assertEqual(len(self.calls('srun')), before)

    def test_collection_requires_comparison_marker_not_only_completed_run(self):
        self.call('submit')
        for task in (1, 2, 3):
            environment = self.task_environment()
            environment['SLURM_ARRAY_TASK_ID'] = str(task)
            self.run_script('slurm/unity/stabilization-validation.sbatch', env=environment)
        next((self.project / 'runs/rostami-01').glob('array-*/results/CPU_GPU_VALIDATION_PASSED')).unlink()
        self.call('collect', expected='failure')

    def test_missing_runtime_and_original_checkout_use_are_rejected(self):
        self.unity_julia.unlink()
        self.call('submit', expected='failure')
        self.call('validate', env={'SATURN_ROSTAMI_RUNTIME_ROOT': self.project / '.unity'}, expected='failure')
        self.assertFalse(self.calls('sbatch'))


if __name__ == '__main__':
    # Select the new subclass explicitly; do not rerun imported baseline tests.
    unittest.main(defaultTest='RostamiTests', verbosity=2)
