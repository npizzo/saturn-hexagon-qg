#!/usr/bin/env python3
"""Preflight one stabilization stage, then reuse the existing immutable launcher.

validate checks provenance and stage hashes without allocating a model,
creating a prepared attempt, or submitting a job. GPU smoke integrations are
the separate validation stage, explicitly submitted with the submit command.
"""

import os
from pathlib import Path
import re
import subprocess
import sys

import experiments_01 as campaign


collect_standard = campaign.collect


REQUIRED_LAUNCH_FILES = {
    'scripts/make_stabilization.jl', 'scripts/make_experiments.jl', 'scripts/run.jl',
    'scripts/validate_stabilization.jl', 'slurm/unity/stabilization-validation.sbatch',
    'slurm/common.sh', 'slurm/ensemble.sbatch', 'slurm/unity/environment.sh',
    'slurm/unity/submit.sh', 'slurm/unity/ensemble.sbatch',
    'slurm/unity/experiments_01.sh', 'slurm/unity/experiments_01.py',
    'slurm/unity/stabilization.sh', 'slurm/unity/stabilization.py',
}


def submit_validation(attempt, expected, project, concurrency, wall_time):
    """Use the inherited frozen receipt format with a separate validation batch.

    The old submission script intentionally has no configurable executable.
    Keep it immutable; only this new validation stage uses this entry point.
    """
    meta = campaign.verify_attempt(attempt, expected, project)
    run_root = Path(os.environ['SATURN_RUN_ROOT']).resolve()
    if run_root != Path(meta['run_root']):
        campaign.fail('run root changed after preparation; restore it or prepare --new-attempt')
    if not Path(os.environ['JULIA_BIN']).is_file() or not os.access(os.environ['JULIA_BIN'], os.X_OK) or not Path(os.environ['JULIA_DEPOT_PATH']).is_dir():
        campaign.fail('Unity Julia environment is missing; validate the existing deployment first')
    environment = {key: value for key, value in os.environ.items() if not key.startswith('SBATCH_')}
    environment.update({'SATURN_EXPERIMENT_ATTEMPT': str(attempt), 'SATURN_EXPERIMENT_INTEGRITY_SHA256': expected,
                        'SATURN_CONFIG_LIST': str(attempt / 'frozen/config-list.tsv'),
                        'SATURN_CONFIG_LIST_SHA256': campaign.digest(attempt / 'frozen/config-list.tsv')})
    (run_root / 'logs').mkdir(parents=True, exist_ok=True)
    log_pattern = str(run_root / 'logs/saturn-ensemble-%A_%a')
    command = ['sbatch', '--parsable', '--export=ALL', '--job-name=saturn-stabilization-validation',
               '--partition=gpu-highvram,uri-gpu', '--constraint=x86_64&a100', '--nodes=1', '--ntasks=1',
               '--gpus=1', '--cpus-per-task=8', '--mem=32G', '--time=' + wall_time, '--no-requeue',
               '--chdir=' + str(project), '--output=' + log_pattern + '.out', '--error=' + log_pattern + '.err',
               '--array=1-{}%{}'.format(meta['count'], concurrency)]
    if os.environ.get('SATURN_UNITY_ACCOUNT'):
        command.append('--account=' + os.environ['SATURN_UNITY_ACCOUNT'])
    command.append(str(project / 'slurm/unity/stabilization-validation.sbatch'))
    campaign.write_new_json(attempt / 'submission-started.json', {
        'command': command, 'started_utc': campaign.datetime.now(campaign.timezone.utc).isoformat(),
        'account_override': os.environ.get('SATURN_UNITY_ACCOUNT'),
        'note': 'Reservation remains even if submission fails or its outcome is unknown.'})
    with (attempt / 'submission.log').open('x') as stream:
        result = subprocess.run(command, env=environment, stdout=stream, stderr=subprocess.STDOUT)
    response = (attempt / 'submission.log').read_text()
    print(response, end='')
    matches = re.findall(r'^([0-9]+)(?:;[^\s;]+)?$', response, re.MULTILINE)
    if result.returncode or len(matches) != 1:
        campaign.fail('submission failed or job ID is uncertain; inspect preserved submission.log before any new attempt: ' + str(attempt))
    campaign.write_new_json(attempt / 'receipt.json', {
        'job_id': matches[0], 'stage': meta['stage'], 'count': meta['count'],
        'concurrency': concurrency, 'time_limit': wall_time, 'integrity_sha256': expected})
    print('Submitted Unity validation job ' + matches[0])
    print('Logs: ' + str(run_root / 'logs'))
    print('Stage receipt: ' + str(attempt / 'receipt.json'))
    print('Only this stage was submitted. Review its results before selecting another stage.')


def collect_validation(attempt):
    # Ordinary run completion can precede the final comparison in the same
    # process. Do not archive that transient state as a passed validation.
    _, _, mapping = campaign.result_map(attempt)
    for row in mapping:
        if len(row['outputs']) != 1:
            campaign.fail('validation needs exactly one completed comparison per member')
        results = Path(row['outputs'][0]['path']) / 'results'
        if not all((results / name).is_file() for name in ('CPU_GPU_VALIDATION_PASSED', 'cpu_gpu_validation.toml')):
            campaign.fail('validation comparison is incomplete or failed: ' + str(results))
    collect_standard(attempt)


def validate_stage(project, suite, stage):
    campaign.verify_provenance(suite / 'provenance.tsv', project)
    listed = {name for _, name in campaign.rows(suite / 'provenance.tsv', 2)}
    missing = REQUIRED_LAUNCH_FILES - listed
    if missing:
        campaign.fail('stabilization provenance is missing: ' + ', '.join(sorted(missing)))
    stages = campaign.rows(suite / 'stages.tsv', 3)
    if len({row[0] for row in stages}) != len(stages):
        campaign.fail('duplicate stage in stages.tsv')
    selected = [row for row in stages if row[0] == stage]
    if len(selected) != 1 or not re.fullmatch('[1-9][0-9]*', selected[0][1]):
        campaign.fail('stage missing or invalid in stages.tsv: ' + stage)
    members = campaign.verify_hashes(suite / stage / 'members.tsv', suite / stage, basenames=True)
    if len(members) != int(selected[0][1]):
        campaign.fail('stage member count does not match stages.tsv')
    # These annotations must be available for the existing immutable preparer.
    campaign.inside(suite, 'campaign.toml')
    campaign.inside(suite, 'physics-audit.tsv')
    print('Validated stabilization stage {}: {} configuration hashes and complete launch provenance.'.format(stage, len(members)), flush=True)


def main(arguments=None):
    arguments = list(sys.argv[1:] if arguments is None else arguments)
    if len(arguments) < 3 or arguments[0] not in ('validate', 'prepare', 'submit', 'status', 'collect'):
        campaign.fail('usage: stabilization.sh validate|prepare|submit|status|collect SUITE STAGE [--time DURATION] [--concurrency N] [--new-attempt]')
    command, suite_arg, stage = arguments[:3]
    if not re.fullmatch('[a-z][a-z0-9_-]*', stage):
        campaign.fail('provide an explicit valid stage name')
    project = Path(os.environ['SATURN_PROJECT']).resolve()
    runtime = Path(os.environ['SATURN_UNITY_ROOT']).resolve()
    if runtime == project / '.unity':
        campaign.fail('use a separate stabilization checkout and reuse the original checkout runtime')
    suite = Path(suite_arg)
    if not suite.is_absolute():
        suite = project / suite
    suite = suite.resolve()
    if command in ('validate', 'prepare', 'submit'):
        validate_stage(project, suite, stage)
    if command == 'validate':
        if len(arguments) != 3:
            campaign.fail('validate accepts only SUITE STAGE; it does not submit the validation stage')
        return
    # status/collect intentionally use frozen receipts even after later source
    # edits. prepare/submit also enforce the existing duplicate/resume guards.
    if command == 'submit' and stage == 'validation':
        campaign.submit = submit_validation
    if command == 'collect' and stage == 'validation':
        campaign.collect = collect_validation
    sys.argv = ['experiments_01.py'] + arguments
    campaign.main()


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError, KeyError) as error:
        print('saturn stabilization: ' + str(error), file=sys.stderr)
        sys.exit(1)
