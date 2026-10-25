#!/usr/bin/env python3
"""Portable, standard-library-only bookkeeping for one Unity experiment stage.

The shell entry point supplies Unity's existing environment. This module never
installs software, chooses a science stage, retries a submission, or cancels jobs.
"""

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile


OPTIONAL_METADATA = ('physics-audit.tsv',)


def fail(message):
    raise ValueError(message)


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def safe_relative(name):
    path = PurePosixPath(name)
    if not name or path.is_absolute() or '..' in path.parts or str(path) != name or '\\' in name:
        fail('unsafe relative path: ' + repr(name))
    return path


def inside(root, relative):
    path = root.joinpath(*safe_relative(relative).parts)
    resolved = path.resolve()
    if root.resolve() not in resolved.parents:
        fail('path escapes its root: ' + str(path))
    if not path.is_file():
        fail('missing file: ' + str(path))
    return path


def rows(path, columns):
    result = []
    for line in path.read_text(encoding='utf-8').splitlines():
        fields = line.split('\t')
        if len(fields) != columns or any(not item or any(ord(c) < 32 for c in item) for item in fields):
            fail('malformed row in ' + str(path))
        result.append(fields)
    if not result:
        fail('empty manifest: ' + str(path))
    return result


def verify_hashes(path, root, basenames=False):
    entries = rows(path, 2)
    seen = set()
    for expected, name in entries:
        safe_relative(name)
        if not re.fullmatch('[0-9a-f]{64}', expected) or name in seen:
            fail('invalid or duplicate hash row in ' + str(path))
        if basenames and ('/' in name or not name.endswith('.toml')):
            fail('member must be a TOML basename: ' + name)
        if digest(inside(root, name)) != expected:
            fail('hash mismatch: ' + str(root / name))
        seen.add(name)
    return entries


def verify_provenance(path, project):
    entries = verify_hashes(path, project)
    names = {name for _, name in entries}
    required = {'Project.toml', 'Manifest.toml', 'scripts/make_experiments.jl'}
    required.update(p.relative_to(project).as_posix() for p in (project / 'src').glob('*.jl'))
    if not required.issubset(names) or not any(name.startswith('src/') for name in names):
        fail('provenance must include all src/*.jl, Project.toml, Manifest.toml, and the generator')


def write_new_json(path, value):
    with path.open('x', encoding='utf-8') as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())


def read_metadata(attempt):
    return json.loads((attempt / 'frozen' / 'attempt.json').read_text())


def verify_frozen(attempt, expected):
    frozen = attempt / 'frozen'
    if not re.fullmatch('[0-9a-f]{64}', expected) or digest(frozen / 'integrity.tsv') != expected:
        fail('frozen integrity manifest changed')
    verify_hashes(frozen / 'integrity.tsv', frozen)
    return read_metadata(attempt)


def verify_attempt(attempt, expected, project):
    frozen = attempt / 'frozen'
    meta = verify_frozen(attempt, expected)
    if meta['project'] != str(project):
        fail('attempt belongs to a different project path; prepare it on this deployment')
    verify_provenance(frozen / 'provenance.tsv', project)
    verify_hashes(frozen / 'members.tsv', frozen / 'configs', basenames=True)
    manifest = rows(frozen / 'config-list.tsv', 2)
    members = rows(frozen / 'members.tsv', 2)
    if manifest != [[h, str(frozen / 'configs' / name)] for h, name in members]:
        fail('absolute configuration list does not match frozen members')
    return meta


def scope_path(suite, stage):
    identity = hashlib.sha256(str(suite).encode('utf-8')).hexdigest()[:20]
    return Path(os.environ['SATURN_UNITY_ROOT']).resolve() / 'experiments-01' / identity / stage


def latest_attempt(scope):
    pointer = scope / 'latest'
    if not pointer.exists():
        return None
    name = pointer.read_text().strip()
    if not re.fullmatch(r'attempt-[A-Za-z0-9_.-]+', name):
        fail('invalid attempt pointer: ' + str(pointer))
    attempt = scope / name
    if not attempt.is_dir():
        fail('missing recorded attempt: ' + str(attempt))
    return attempt


def prepare(suite, stage, project, new_attempt):
    scope = scope_path(suite, stage)
    prior = latest_attempt(scope)
    if prior and not new_attempt:
        if (prior / 'submission-started.json').exists():
            fail('stage already submitted or submission outcome uncertain; inspect status before using --new-attempt: ' + str(prior))
        expected = (prior / 'integrity.sha256').read_text().strip()
        verify_attempt(prior, expected, project)
        # Never silently use an older frozen stage after the suite was edited.
        for name in ('members.tsv',):
            if digest(suite / stage / name) != digest(prior / 'frozen' / name):
                fail('suite stage changed after preparation; use --new-attempt')
        for name in ('provenance.tsv', 'campaign.toml', 'stages.tsv'):
            if digest(suite / name) != digest(prior / 'frozen' / name):
                fail('suite metadata changed after preparation; use --new-attempt')
        for name in OPTIONAL_METADATA:
            original, frozen_copy = suite / name, prior / 'frozen' / name
            if original.exists() != frozen_copy.exists() or (original.exists() and digest(original) != digest(frozen_copy)):
                fail('optional suite metadata changed after preparation; use --new-attempt: ' + name)
        verify_hashes(suite / stage / 'members.tsv', suite / stage, basenames=True)
        return prior, expected
    stages = rows(suite / 'stages.tsv', 3)
    if len({row[0] for row in stages}) != len(stages):
        fail('duplicate stage in stages.tsv')
    selected = [row for row in stages if row[0] == stage]
    if len(selected) != 1 or not re.fullmatch('[1-9][0-9]*', selected[0][1]):
        fail('stage missing or invalid in stages.tsv: ' + stage)
    verify_provenance(suite / 'provenance.tsv', project)
    members = verify_hashes(suite / stage / 'members.tsv', suite / stage, basenames=True)
    if len(members) != int(selected[0][1]):
        fail('stage member count does not match stages.tsv')
    scope.mkdir(parents=True, exist_ok=True)
    if new_attempt:
        attempt = Path(tempfile.mkdtemp(prefix='attempt-' + datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S') + '-', dir=str(scope)))
    else:
        attempt = scope / 'attempt-001'
        attempt.mkdir()  # Exclusive reservation also handles concurrent prepares.
    frozen = attempt / 'frozen'
    (frozen / 'configs').mkdir(parents=True)
    for _, name in members:
        shutil.copyfile(suite / stage / name, frozen / 'configs' / name)
    shutil.copyfile(suite / stage / 'members.tsv', frozen / 'members.tsv')
    for name in ('provenance.tsv', 'campaign.toml', 'stages.tsv'):
        shutil.copyfile(suite / name, frozen / name)
    for name in OPTIONAL_METADATA:
        if (suite / name).exists():
            shutil.copyfile(inside(suite, name), frozen / name)
    metadata = {'schema': 1, 'suite': str(suite), 'stage': stage, 'count': len(members),
                'project': str(project), 'run_root': str(Path(os.environ['SATURN_RUN_ROOT']).resolve()),
                'review_note': selected[0][2], 'prepared_utc': datetime.now(timezone.utc).isoformat()}
    write_new_json(frozen / 'attempt.json', metadata)
    (frozen / 'config-list.tsv').write_text(''.join(h + '\t' + str(frozen / 'configs' / name) + '\n' for h, name in members))
    files = sorted(p for p in frozen.rglob('*') if p.is_file())
    (frozen / 'integrity.tsv').write_text(''.join(digest(p) + '\t' + p.relative_to(frozen).as_posix() + '\n' for p in files))
    expected = digest(frozen / 'integrity.tsv')
    (attempt / 'integrity.sha256').write_text(expected + '\n')
    verify_attempt(attempt, expected, project)  # Detect changes during copying too.
    for path in frozen.rglob('*'):
        path.chmod(0o555 if path.is_dir() else 0o444)
    frozen.chmod(0o555)
    with tempfile.NamedTemporaryFile(mode='w', dir=str(scope), delete=False) as stream:
        stream.write(attempt.name + '\n')
        pointer = Path(stream.name)
    os.replace(str(pointer), str(scope / 'latest'))
    return attempt, expected


def submit(attempt, expected, project, concurrency, wall_time):
    meta = verify_attempt(attempt, expected, project)
    if Path(os.environ['SATURN_RUN_ROOT']).resolve() != Path(meta['run_root']):
        fail('run root changed after preparation; restore it or prepare --new-attempt')
    if not Path(os.environ['JULIA_BIN']).is_file() or not os.access(os.environ['JULIA_BIN'], os.X_OK) or not Path(os.environ['JULIA_DEPOT_PATH']).is_dir():
        fail('Unity Julia environment is missing; validate the existing deployment first')
    environment = dict(os.environ)
    environment.update({'SATURN_EXPERIMENT_ATTEMPT': str(attempt), 'SATURN_EXPERIMENT_INTEGRITY_SHA256': expected,
                        'SATURN_CONFIG_LIST': str(attempt / 'frozen' / 'config-list.tsv'),
                        'SATURN_CONFIG_LIST_SHA256': digest(attempt / 'frozen' / 'config-list.tsv')})
    command = ['bash', str(project / 'slurm/unity/submit.sh'), 'ensemble',
               '--array=1-{}%{}'.format(meta['count'], concurrency), '--time=' + wall_time]
    write_new_json(attempt / 'submission-started.json', {'command': command, 'started_utc': datetime.now(timezone.utc).isoformat(),
                                                       'account_override': os.environ.get('SATURN_UNITY_ACCOUNT'),
                                                       'note': 'Reservation remains even if submission fails or its outcome is unknown.'})
    # Preserve output even if interrupted after sbatch accepted the job. Never
    # infer this stage's job from the shared, replaceable last-job-id pointer.
    with (attempt / 'submission.log').open('x') as stream:
        result = subprocess.run(command, env=environment, stdout=stream, stderr=subprocess.STDOUT)
    response = (attempt / 'submission.log').read_text()
    print(response, end='')
    matches = re.findall(r'^Submitted Unity ensemble job ([0-9]+)$', response, re.MULTILINE)
    if result.returncode or len(matches) != 1:
        fail('submission failed or job ID is uncertain; inspect preserved submission.log before any new attempt: ' + str(attempt))
    write_new_json(attempt / 'receipt.json', {'job_id': matches[0], 'stage': meta['stage'], 'count': meta['count'],
                                           'concurrency': concurrency, 'time_limit': wall_time, 'integrity_sha256': expected})
    print('Stage receipt: ' + str(attempt / 'receipt.json'))
    print('Only this stage was submitted. Review its results before selecting another stage.')


def result_map(attempt):
    receipt = json.loads((attempt / 'receipt.json').read_text())
    meta = verify_frozen(attempt, receipt['integrity_sha256'])
    if receipt['stage'] != meta['stage'] or receipt['count'] != meta['count']:
        fail('receipt does not match frozen stage metadata')
    job = receipt['job_id']
    if not re.fullmatch('[0-9]+', job):
        fail('invalid job ID in receipt')
    run_root = Path(meta['run_root'])
    mapping = []
    for index, (expected, name) in enumerate(rows(attempt / 'frozen' / 'members.tsv', 2), 1):
        paths = sorted(run_root.glob('array-{}-{}.*'.format(job, index)))
        mapping.append({'task': index, 'config': name, 'config_sha256': expected,
                        'outputs': [{'path': str(path), 'completed': (path / 'results/COMPLETED').is_file(),
                                     'config_matches': (path / 'config.toml').is_file() and digest(path / 'config.toml') == expected} for path in paths],
                        'logs': [str(run_root / 'logs' / 'saturn-ensemble-{}_{}.{}'.format(job, index, extension)) for extension in ('out', 'err')]})
    return meta, receipt, mapping


def status(attempt):
    print('Attempt: ' + str(attempt))
    if not (attempt / 'receipt.json').is_file():
        print('No confirmed job receipt. ' + ('Inspect submission.log before retrying.' if (attempt / 'submission-started.json').exists() else 'Prepared; not submitted.'))
        return
    _, receipt, mapping = result_map(attempt)
    print('Stage {}: job {}; {} members'.format(receipt['stage'], receipt['job_id'], receipt['count']))
    for row in mapping:
        print('{}: {} | {}'.format(row['task'], row['config'], ', '.join(p['path'] + (' [COMPLETED]' if p['completed'] and p['config_matches'] else ' [incomplete or config mismatch]') for p in row['outputs']) or 'no output directory yet'))
    for command in (['squeue', '-j', receipt['job_id'], '-o', '%.18i %.30j %.12T %.10M %R'],
                    ['sacct', '-j', receipt['job_id'], '--format=JobID,JobName,State,ExitCode,Elapsed,Reason']):
        result = subprocess.run(command)
        if result.returncode:
            print('{} returned {}; inspect the other scheduler output too.'.format(command[0], result.returncode), file=sys.stderr)


def collect(attempt):
    meta, receipt, mapping = result_map(attempt)
    target = attempt / ('results-' + receipt['job_id'] + '.tar.gz')
    if target.exists():
        fail('collection already exists; refusing overwrite: ' + str(target))
    # Completion and configuration identity come from each recorded output, not
    # scheduler state or directory names. They do not establish scientific validity.
    collection = {'receipt': receipt, 'members': mapping, 'collected_utc': datetime.now(timezone.utc).isoformat(),
                  'all_members_completed': all(len(r['outputs']) == 1 and r['outputs'][0]['completed'] and r['outputs'][0]['config_matches'] for r in mapping)}
    if not collection['all_members_completed']:
        fail('some tasks lack one completed output with the expected configuration; inspect status before collection')
    fd, temporary = tempfile.mkstemp(prefix='.collection-', dir=str(attempt))
    os.close(fd)
    try:
        with tarfile.open(temporary, 'w:gz') as archive:
            archive.add(str(attempt / 'frozen'), arcname='inputs', recursive=True)
            for name in ('integrity.sha256', 'submission-started.json', 'submission.log', 'receipt.json'):
                archive.add(str(attempt / name), arcname='submission/' + name)
            encoded = (json.dumps(collection, indent=2, sort_keys=True) + '\n').encode('utf-8')
            import io
            info = tarfile.TarInfo('collection.json')
            info.size = len(encoded)
            archive.addfile(info, io.BytesIO(encoded))
            for row in mapping:
                for output in row['outputs']:
                    path = Path(output['path'])
                    if path.is_symlink():
                        fail('refusing a symlink output directory: ' + str(path))
                    archive.add(str(path), arcname='runs/' + path.name)
                for filename in row['logs']:
                    path = Path(filename)
                    if path.is_file():
                        archive.add(str(path), arcname='logs/' + path.name)
        os.link(temporary, str(target))  # Atomic publication without overwrite.
    finally:
        os.unlink(temporary)
    print('Collected completed stage: ' + str(target))
    print('SHA256: ' + digest(target))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('prepare', 'submit', 'status', 'collect', 'verify-task'))
    parser.add_argument('suite', nargs='?')
    parser.add_argument('stage', nargs='?')
    parser.add_argument('--new-attempt', action='store_true', help='explicitly reserve another attempt; never cancels prior jobs')
    parser.add_argument('--concurrency', type=int, default=2)
    parser.add_argument('--time', default='01:00:00', help='Slurm wall-time limit; default one hour')
    args = parser.parse_args()
    project = Path(os.environ['SATURN_PROJECT']).resolve()
    if args.command == 'verify-task':
        attempt = Path(os.environ['SATURN_EXPERIMENT_ATTEMPT'])
        if not attempt.is_absolute():
            fail('experiment attempt path must be absolute')
        verify_attempt(attempt, os.environ['SATURN_EXPERIMENT_INTEGRITY_SHA256'], project)
        if os.environ.get('SATURN_CONFIG_LIST') != str(attempt / 'frozen/config-list.tsv'):
            fail('array config list does not belong to this experiment attempt')
        print('Experiment source, dependency, and configuration hashes verified.')
        return
    if args.suite is None or args.stage is None or not re.fullmatch('[a-z][a-z0-9_-]*', args.stage):
        parser.error('provide SUITE and an explicit STAGE')
    if args.concurrency < 1 or not re.fullmatch(r'[0-9]+(?:[-:][0-9]+){0,3}', args.time):
        parser.error('concurrency must be positive and time must be a numeric Slurm duration')
    if args.command in ('prepare', 'submit') and (os.environ.get('SLURM_JOB_ID') or os.environ.get('SATURN_RESUME')):
        fail('prepare/submit from a login shell outside an allocation, with SATURN_RESUME unset')
    suite = Path(args.suite)
    if not suite.is_absolute():
        suite = project / suite
    suite = suite.resolve()
    if args.command in ('prepare', 'submit'):
        attempt, expected = prepare(suite, args.stage, project, args.new_attempt)
        print('Prepared attempt: ' + str(attempt), flush=True)
        if args.command == 'submit':
            submit(attempt, expected, project, args.concurrency, args.time)
    else:
        if args.new_attempt:
            parser.error('--new-attempt is valid only with prepare or submit')
        attempt = latest_attempt(scope_path(suite, args.stage))
        if attempt is None:
            fail('no prepared attempt for this suite and stage')
        if args.command == 'status':
            status(attempt)
        else:
            collect(attempt)


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError, KeyError) as error:
        print('saturn experiments: ' + str(error), file=sys.stderr)
        sys.exit(1)
