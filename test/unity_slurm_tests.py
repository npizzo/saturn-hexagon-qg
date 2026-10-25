#!/usr/bin/env python3
"""Offline integration checks for the Unity deployment entry points.

Run with ``python3 test/unity_slurm_tests.py``. Every subprocess runs in a
temporary project whose path contains spaces. Scheduler, runtime, and download
commands are mocks; these tests never submit a cluster job or use the network.
The mocks retain argument boundaries and selected environment variables, while
real shell scripts perform their own path handling, staging, and hash checks.
"""

import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest


PROJECT = Path(__file__).resolve().parents[1]

MOCK_TOOL = r'''
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

name = Path(sys.argv[0]).name
args = sys.argv[1:]
entry = {"tool": name, "args": args, "env": {
    k: v for k, v in os.environ.items()
    if k.startswith(("SATURN_", "SBATCH_", "SLURM_", "JULIA_", "OPENBLAS_"))
}}
with open(os.environ["MOCK_LOG"], "a", encoding="utf-8") as stream:
    stream.write(json.dumps(entry) + "\n")
if name == os.environ.get("MOCK_FAIL_TOOL") and os.environ.get("MOCK_FAIL_MATCH", "") in " ".join(args):
    print("intentional mock failure: " + name, file=sys.stderr)
    sys.exit(37)
if name == "uname":
    print("Linux" if args == ["-s"] else os.environ.get("MOCK_ARCH", "x86_64"))
elif name == "sbatch":
    if os.environ.get("MOCK_SBATCH_WARNING") == "1":
        print("sbatch: harmless scheduler notice", file=sys.stderr)
    number = 12345
    if os.environ.get("MOCK_SBATCH_SEQUENCE") == "1":
        records = [json.loads(line) for line in Path(os.environ["MOCK_LOG"]).read_text().splitlines()]
        number += sum(r["tool"] == "sbatch" for r in records) - 1
    print(str(number) if "--parsable" in args else "Submitted batch job " + str(number))
elif name in ("squeue", "sacct"):
    print("offline scheduler status")
elif name == "srun":
    command = next(i for i, value in enumerate(args) if not value.startswith("-"))
    sys.exit(subprocess.run(args[command:]).returncode)
elif name == "julia":
    if "--version" in args:
        print("julia version 1.10.12")
elif name == "curl":
    output = args[args.index("--output") + 1]
    Path(output).write_bytes(b"offline mock archive; never a real download")
elif name == "tar":
    stage = Path(args[args.index("-C") + 1])
    executable = stage / "julia-1.10.12" / "bin" / "julia"
    executable.parent.mkdir(parents=True)
    shutil.copy2(os.environ["MOCK_JULIA"], executable)
elif name == "sha256sum":
    if "--check" in args:
        # Archive checks cannot verify mock bytes against a production digest;
        # record the exact requested checksum and test the failure branch.
        line = sys.stdin.read()
        with open(os.environ["MOCK_LOG"], "a", encoding="utf-8") as stream:
            stream.write(json.dumps({"tool": "checksum-input", "text": line}) + "\n")
        if os.environ.get("MOCK_BAD_CHECKSUM") == "1":
            print("mock checksum mismatch", file=sys.stderr)
            sys.exit(1)
    else:
        for filename in args:
            print(hashlib.sha256(Path(filename).read_bytes()).hexdigest() + "  " + filename)
else:
    raise RuntimeError("unexpected mock executable " + name)
'''


class SlurmHarness(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="saturn Unity tests ")
        self.addCleanup(self.temporary.cleanup)
        # Match shell `pwd -P` (macOS /var is a symlink to /private/var).
        self.root = Path(self.temporary.name).resolve()
        self.project = self.root / "project with spaces"
        self.project.mkdir()
        shutil.copytree(PROJECT / "slurm", self.project / "slurm")
        for name in ("Project.toml", "Manifest.toml"):
            (self.project / name).write_text("# offline shell-test fixture\n")
        for name in ("scripts/run.jl", "scripts/validate_gpu.jl", "scripts/validate_gpu_restart.jl", "test/runtests.jl"):
            destination = self.project / name
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_text("# mocked Julia invocation\n")
        (self.project / "configs").mkdir()
        self.config = self.project / "configs" / "smoke.toml"
        self.config.write_text('[experiment]\nname = "offline fixture"\n')
        self.spaced_config = self.project / "configs" / "pilot with spaces.toml"
        self.spaced_config.write_text('[experiment]\nname = "spaced fixture"\n')
        (self.project / "configs/unity_benchmark.toml").write_text('[experiment]\nname = "benchmark fixture"\n')
        followup = self.project / "configs/unity_followup_01"
        followup.mkdir()
        for name in ("maintenance_forced", "maintenance_unforced"):
            (followup / (name + ".toml")).write_text('[experiment]\nname = "' + name + '"\n')
        self.bin = self.root / "mock tools"
        self.bin.mkdir()
        for name in ("sbatch", "squeue", "sacct", "srun", "julia", "curl", "tar", "sha256sum", "uname"):
            executable = self.bin / name
            executable.write_text("#!" + sys.executable + "\n" + MOCK_TOOL)
            executable.chmod(0o755)
        self.log = self.root / "tool calls.jsonl"
        self.run_root = self.root / "run output"
        self.depot = self.root / "Julia depot"
        self.depot.mkdir()
        self.runtime = self.root / "Julia runtime"
        self.unity_root = self.project / ".unity"
        self.unity_julia = self.unity_root / "runtime" / "julia-1.10.12" / "bin" / "julia"
        self.unity_julia.parent.mkdir(parents=True)
        shutil.copy2(self.bin / "julia", self.unity_julia)
        (self.unity_root / "julia-depot").mkdir()
        self.env = {k: v for k, v in os.environ.items()
                    if not k.startswith(("SATURN_", "SBATCH_", "SLURM_", "JULIA_", "MOCK_"))}
        self.env.update({
            "PATH": str(self.bin) + os.pathsep + os.environ.get("PATH", "/usr/bin:/bin"),
            "MOCK_LOG": str(self.log),
            "MOCK_JULIA": str(self.bin / "julia"),
            "SATURN_PROJECT": str(self.project),
            "SATURN_RUN_ROOT": str(self.run_root),
            "JULIA_BIN": str(self.bin / "julia"),
            "JULIA_DEPOT_PATH": str(self.depot),
            "SLURM_JOB_ID": "12345",
            "SLURM_JOB_ACCOUNT": "unity-default",
            "SLURM_CPUS_PER_TASK": "8",
        })

    def run_script(self, relative, *args, env=None, expected=0):
        environment = dict(self.env)
        if env:
            for key, value in env.items():
                if value is None:
                    environment.pop(key, None)
                else:
                    environment[key] = str(value)
        result = subprocess.run(
            ["bash", str(self.project / relative), *map(str, args)],
            cwd=self.root, env=environment, text=True, capture_output=True,
            # Several mocked processes may start under a busy desktop host.
            # This bounds a hung check without mistaking startup delay for failure.
            timeout=60, check=False,
        )
        if expected == 0:
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        elif expected == "failure":
            self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        else:
            self.assertEqual(result.returncode, expected, result.stdout + result.stderr)
        return result

    def calls(self, tool=None):
        if not self.log.exists():
            return []
        records = [json.loads(line) for line in self.log.read_text().splitlines()]
        return records if tool is None else [r for r in records if r["tool"] == tool]

    def reset_calls(self):
        self.log.unlink(missing_ok=True)

    def manifest(self):
        path = self.root / "frozen manifests" / "config list.tsv"
        self.run_script("slurm/prepare_array.sh", path, self.config, self.spaced_config)
        self.reset_calls()
        return path, hashlib.sha256(path.read_bytes()).hexdigest()

    def array_environment(self, manifest, digest, task="2"):
        return {"SATURN_CONFIG_LIST": manifest,
                "SATURN_CONFIG_LIST_SHA256": digest,
                "SLURM_ARRAY_JOB_ID": "12000",
                "SLURM_ARRAY_TASK_ID": task}

    def assert_gpu_launches(self, count=1):
        launches = [r for r in self.calls("julia")
                    if str(self.project / "scripts/run.jl") in r["args"]]
        self.assertEqual(len(launches), count, self.calls())
        for launch in launches:
            index = launch["args"].index("--device")
            self.assertEqual(launch["args"][index + 1], "gpu")
        return launches


class UnitySubmissionTests(SlurmHarness):
    def setUp(self):
        super().setUp()
        # Submission is a login-node action, outside an existing allocation.
        self.env.pop("SLURM_JOB_ID")
        self.env.pop("SLURM_JOB_ACCOUNT")
        self.env.pop("SLURM_CPUS_PER_TASK")

    def test_default_submission_discards_inherited_delta_settings(self):
        self.run_script("slurm/unity/submit.sh", "validate", env={
            "SBATCH_ACCOUNT": "bhcr",
            "SBATCH_QOS": "short",
            "SBATCH_ARRAY_INX": "9-10",
            "SATURN_PROJECT": "/old/delta/project",
            "SATURN_RUN_ROOT": "/old/delta/output",
            "JULIA_BIN": "/old/delta/aarch64/julia",
            "JULIA_DEPOT_PATH": "/old/delta/depot",
        })
        calls = self.calls("sbatch")
        self.assertEqual(len(calls), 1)
        args, env = calls[0]["args"], calls[0]["env"]
        self.assertNotIn("SBATCH_ACCOUNT", env)
        self.assertNotIn("SBATCH_QOS", env)
        self.assertNotIn("SBATCH_ARRAY_INX", env)
        self.assertFalse(any(key.startswith("SBATCH_") for key in env))
        self.assertNotIn("bhcr", " ".join(args))
        self.assertFalse(any(a == "-A" or a.startswith("--account") for a in args))
        self.assertIn("--partition=gpu-highvram,uri-gpu", args)
        self.assertIn("--constraint=x86_64&a100", args)
        self.assertIn("--gpus=1", args)
        self.assertIn("--no-requeue", args)
        self.assertFalse(any(a.startswith("--qos") for a in args))
        self.assertEqual(env["SATURN_PROJECT"], str(self.project))
        self.assertEqual(env["SATURN_RUN_ROOT"], str(self.project / "runs" / "unity"))
        self.assertEqual(env["JULIA_BIN"], str(self.unity_julia))
        self.assertEqual(env["JULIA_DEPOT_PATH"], str(self.unity_root / "julia-depot"))
        self.assertEqual(env["JULIA_CPU_TARGET"], "generic")
        self.assertEqual(env["JULIA_NUM_THREADS"], "8")
        self.assertEqual(env["JULIA_NUM_PRECOMPILE_TASKS"], "4")
        self.assertEqual(env["OPENBLAS_NUM_THREADS"], "1")
        self.assertEqual((self.unity_root / "last-job-id").read_text().strip(), "12345")
        self.assertFalse(self.calls("srun"))

    def test_run_submission_preserves_config_and_explicit_overrides(self):
        options = ["--account=cli-account", "--partition=chosen-gpu", "--time=00:07:00", "--qos=chosen-qos",
                   "--cpus-per-task=4", "--mem=12G", "--constraint=x86_64&a100"]
        custom_root = self.root / "custom Unity runtime root"
        custom_runs = self.root / "custom run destination"
        custom_julia = custom_root / "runtime/julia-1.10.12/bin/julia"
        custom_julia.parent.mkdir(parents=True)
        shutil.copy2(self.bin / "julia", custom_julia)
        (custom_root / "julia-depot").mkdir()
        self.run_script("slurm/unity/submit.sh", "run", self.spaced_config, *options, env={
            "SATURN_UNITY_ROOT": custom_root,
            "SATURN_UNITY_RUN_ROOT": custom_runs,
            "SATURN_UNITY_ACCOUNT": "environment-account",
            "SBATCH_ACCOUNT": "bhcr",
            "SBATCH_QOS": "short",
            "SBATCH_ARRAY_INX": "9-10",
        })
        call = self.calls("sbatch")[0]
        args, env = call["args"], call["env"]
        self.assertEqual(args.count(str(self.spaced_config)), 1)
        script = str(self.project / "slurm/unity/single_gpu.sbatch")
        self.assertIn(script, args)
        self.assertGreater(args.index(str(self.spaced_config)), args.index(script))
        for option in options:
            self.assertIn(option, args)
            self.assertLess(args.index(option), args.index(script))
        self.assertIn("--account=environment-account", args)
        self.assertGreater(args.index("--account=cli-account"), args.index("--account=environment-account"))
        self.assertNotIn("SBATCH_ACCOUNT", env)
        self.assertNotIn("SBATCH_QOS", env)
        self.assertNotIn("SBATCH_ARRAY_INX", env)
        self.assertEqual(env["SATURN_RUN_ROOT"], str(custom_runs))
        self.assertEqual(env["JULIA_BIN"], str(custom_root / "runtime/julia-1.10.12/bin/julia"))
        self.assertEqual(env["JULIA_DEPOT_PATH"], str(custom_root / "julia-depot"))
        self.assertEqual((custom_root / "last-job-id").read_text().strip(), "12345")

    def test_ensemble_submission_keeps_manifest_and_array_limit(self):
        manifest, digest = self.manifest()
        self.run_script("slurm/unity/submit.sh", "ensemble", "--array=1-2%1",
                        env=self.array_environment(manifest, digest))
        call = self.calls("sbatch")[0]
        self.assertIn("--array=1-2%1", call["args"])
        self.assertIn(str(self.project / "slurm/unity/ensemble.sbatch"), call["args"])
        self.assertEqual(call["env"]["SATURN_CONFIG_LIST"], str(manifest))
        self.assertEqual(call["env"]["SATURN_CONFIG_LIST_SHA256"], digest)

    def test_failed_submission_does_not_record_a_job(self):
        self.run_script("slurm/unity/submit.sh", "validate",
                        env={"MOCK_FAIL_TOOL": "sbatch"}, expected="failure")
        self.assertFalse((self.unity_root / "last-job-id").exists())
        self.assertFalse(self.calls("srun"))

    def test_non_x86_unity_host_is_rejected_before_submission(self):
        self.run_script("slurm/unity/submit.sh", "validate",
                        env={"MOCK_ARCH": "aarch64"}, expected="failure")
        self.assertFalse(self.calls("sbatch"))
        self.assertFalse(self.calls("curl"))

    def test_uninstalled_runtime_is_not_submitted_for_production(self):
        self.unity_julia.unlink()
        self.run_script("slurm/unity/submit.sh", "run", self.spaced_config, expected="failure")
        self.assertFalse(self.calls("sbatch"))
        self.assertFalse(self.calls("curl"))

    def test_restart_submission_and_afterok_dependency_preserve_arguments(self):
        self.run_script("slurm/unity/submit.sh", "restart-check", self.spaced_config)
        call = self.calls("sbatch")[0]
        script = str(self.project / "slurm/unity/validate_gpu_restart.sbatch")
        self.assertEqual(call["args"][-2:], [script, str(self.spaced_config)])
        self.assertIn("--job-name=saturn-restart-check", call["args"])
        self.assertFalse(self.calls("srun"))
        self.reset_calls()
        self.run_script("slurm/unity/submit.sh", "run", self.spaced_config, "--dependency=afterok:12345")
        args = self.calls("sbatch")[0]["args"]
        self.assertIn("--dependency=afterok:12345", args)
        self.assertLess(args.index("--dependency=afterok:12345"), args.index(str(self.project / "slurm/unity/single_gpu.sbatch")))

    def test_restart_submission_rejects_existing_resume_and_missing_runtime(self):
        self.run_script("slurm/unity/submit.sh", "restart-check", self.config,
                        env={"SATURN_RESUME": "/unrelated/checkpoint.jld2"}, expected="failure")
        self.assertFalse(self.calls("sbatch"))
        self.unity_julia.unlink()
        self.run_script("slurm/unity/submit.sh", "restart-check", self.config, expected="failure")
        self.assertFalse(self.calls("sbatch"))

    def test_followup_records_three_jobs_and_gates_both_maintenance_runs(self):
        self.run_script("slurm/unity/followup_01.sh", "submit", env={"MOCK_SBATCH_SEQUENCE": "1", "MOCK_SBATCH_WARNING": "1"})
        calls = self.calls("sbatch")
        self.assertEqual(len(calls), 3)
        self.assertFalse(any(a.startswith("--dependency") for a in calls[0]["args"]))
        for call in calls[1:]:
            self.assertIn("--dependency=afterok:12345", call["args"])
        ledger = self.unity_root / "followup-01-jobs.tsv"
        expected = "role\tjob_id\nrestart_check\t12345\nmaintenance_forced\t12346\nmaintenance_unforced\t12347\n"
        self.assertEqual(ledger.read_text(), expected)
        self.reset_calls()
        self.run_script("slurm/unity/followup_01.sh", "status")
        self.assertEqual(self.calls("squeue")[0]["args"][1], "12345,12346,12347")
        self.assertEqual(self.calls("sacct")[0]["args"][1], "12345,12346,12347")
        self.assertFalse(self.calls("sbatch"))
        self.assertEqual(ledger.read_text(), expected)
        self.reset_calls()
        self.run_script("slurm/unity/followup_01.sh", "submit", expected="failure")
        self.assertFalse(self.calls("sbatch"))
        self.assertEqual(ledger.read_text(), expected)

    def test_followup_preserves_partial_submission_and_status_uses_only_saved_ids(self):
        self.run_script("slurm/unity/followup_01.sh", "submit", env={
            "MOCK_FAIL_TOOL": "sbatch", "MOCK_FAIL_MATCH": "maintenance_forced.toml",
            "MOCK_SBATCH_SEQUENCE": "1",
        }, expected="failure")
        self.assertEqual(len(self.calls("sbatch")), 2)
        ledger = self.unity_root / "followup-01-jobs.tsv"
        self.assertEqual(ledger.read_text(), "role\tjob_id\nrestart_check\t12345\n")
        self.assertIn("intentional mock failure", (self.unity_root / "followup-01-submission.log").read_text())
        # The generic last-job pointer is deliberately unrelated to this ledger.
        (self.unity_root / "last-job-id").write_text("99999\n")
        self.reset_calls()
        self.run_script("slurm/unity/followup_01.sh", "status")
        self.assertEqual(self.calls("squeue")[0]["args"][1], "12345")
        self.assertEqual(self.calls("sacct")[0]["args"][1], "12345")
        self.run_script("slurm/unity/followup_01.sh", "submit", expected="failure")
        self.assertFalse(self.calls("sbatch"))

    def test_followup_status_still_queries_accounting_when_queue_is_empty(self):
        ledger = self.unity_root / "followup-01-jobs.tsv"
        ledger.write_text("role\tjob_id\nrestart_check\t12345\n")
        self.run_script("slurm/unity/followup_01.sh", "status", env={"MOCK_FAIL_TOOL": "squeue"})
        self.assertEqual(len(self.calls("sacct")), 1)
        self.assertFalse(self.calls("sbatch"))


class UnityExecutionTests(SlurmHarness):
    def test_restart_check_reuses_installation_and_stages_config(self):
        self.run_script("slurm/unity/validate_gpu_restart.sbatch", self.spaced_config)
        launches = [r for r in self.calls("julia")
                    if str(self.project / "scripts/validate_gpu_restart.jl") in r["args"]]
        self.assertEqual(len(launches), 1)
        args = launches[0]["args"]
        self.assertEqual(args[args.index("--device") + 1], "gpu")
        staged = Path(args[args.index(str(self.project / "scripts/validate_gpu_restart.jl")) + 1])
        self.assertEqual(staged.read_bytes(), self.spaced_config.read_bytes())
        self.assertEqual(args[args.index("--output") + 1], str(staged.parent / "results"))
        self.assertTrue((staged.parent / "job-context.txt").is_file())
        self.assertEqual(len(self.calls("srun")), 1)
        self.assertFalse(self.calls("curl"))
        self.assertFalse(any(any("Pkg.instantiate()" in a for a in r["args"]) for r in self.calls("julia")))

    def test_restart_failure_propagates_without_other_work(self):
        self.run_script("slurm/unity/validate_gpu_restart.sbatch", self.config, env={
            "MOCK_FAIL_TOOL": "julia", "MOCK_FAIL_MATCH": "scripts/validate_gpu_restart.jl",
        }, expected=37)
        self.assertEqual(len(self.calls("srun")), 1)
        self.assert_gpu_launches(0)
        self.assertFalse(self.calls("curl"))

    def test_restart_execution_rejects_resume_before_launch(self):
        self.run_script("slurm/unity/validate_gpu_restart.sbatch", self.config, env={
            "SATURN_RESUME": "/unrelated/checkpoint.jld2",
        }, expected="failure")
        self.assertFalse(self.calls("srun"))
        self.assertFalse(self.calls("julia"))

    def test_validation_bootstraps_once_then_checks_and_runs_gpu(self):
        self.run_script("slurm/unity/validate_gpu.sbatch")
        bootstraps = [r for r in self.calls("srun")
                      if str(self.project / "slurm/bootstrap.sh") in r["args"]]
        self.assertEqual(len(bootstraps), 1)
        installs = [r for r in self.calls("julia") if any("Pkg.instantiate()" in a for a in r["args"])]
        self.assertEqual(len(installs), 1)
        self.assertEqual(installs[0]["env"]["JULIA_DEPOT_PATH"], str(self.unity_root / "julia-depot"))
        self.assertEqual(sum(str(self.project / "test/runtests.jl") in r["args"] for r in self.calls("julia")), 1)
        self.assertEqual(sum(str(self.project / "scripts/validate_gpu.jl") in r["args"] for r in self.calls("julia")), 1)
        self.assert_gpu_launches()
        self.assertFalse(self.calls("curl"))  # Existing pinned runtime is reused.

    def test_validation_bootstrap_failure_stops_all_later_phases(self):
        self.run_script("slurm/unity/validate_gpu.sbatch", env={
            "MOCK_FAIL_TOOL": "julia", "MOCK_FAIL_MATCH": "Pkg.instantiate()",
        }, expected="failure")
        self.assertEqual(len(self.calls("srun")), 1)
        self.assertFalse(any(str(self.project / "test/runtests.jl") in r["args"] for r in self.calls("julia")))
        self.assert_gpu_launches(0)

    def test_validation_failure_does_not_launch_production_or_cpu_fallback(self):
        self.run_script("slurm/unity/validate_gpu.sbatch", env={
            "MOCK_FAIL_TOOL": "julia", "MOCK_FAIL_MATCH": "scripts/validate_gpu.jl",
        }, expected="failure")
        self.assert_gpu_launches(0)
        self.assertFalse(any("cpu" in r["args"] for r in self.calls("julia")))

    def test_validation_resume_is_rejected_before_bootstrap(self):
        self.run_script("slurm/unity/validate_gpu.sbatch", env={
            "SATURN_RESUME": str(self.root / "unrelated checkpoint.jld2"),
        }, expected="failure")
        self.assertFalse(self.calls("srun"))
        self.assertFalse(self.calls("julia"))
        self.assertFalse(self.calls("curl"))

    def test_single_run_never_bootstraps_and_does_not_fallback(self):
        self.run_script("slurm/unity/single_gpu.sbatch", self.spaced_config, env={
            "MOCK_FAIL_TOOL": "julia", "MOCK_FAIL_MATCH": "scripts/run.jl",
            "SLURM_CPUS_PER_TASK": "4", "JULIA_NUM_THREADS": "64",
        }, expected="failure")
        launch = self.assert_gpu_launches()[0]
        self.assertEqual(launch["env"]["JULIA_NUM_THREADS"], "4")
        self.assertEqual(len(self.calls("srun")), 1)
        self.assertFalse(self.calls("curl"))
        self.assertFalse(any(any("Pkg.instantiate()" in a for a in r["args"]) for r in self.calls("julia")))
        self.assertFalse(any("cpu" in r["args"] for r in self.calls("julia")))


class BootstrapTests(SlurmHarness):
    def test_pinned_architecture_download_selection_and_reuse(self):
        expected_hashes = {
            "x86_64": "b0c50dfe349d76ffc8580d24e0c0be623988a6d91ed57ad8b73c0137fd27bc0d",
            "aarch64": "cc257b21967049642bfe17de1aff7dae17477eba4f423f98cb400c002080e5a2",
        }
        for arch in ("x86_64", "aarch64"):
            with self.subTest(arch=arch):
                self.reset_calls()
                runtime = self.runtime / arch
                self.run_script("slurm/bootstrap.sh", runtime, self.depot, env={"MOCK_ARCH": arch})
                downloads = self.calls("curl")
                self.assertEqual(len(downloads), 1)
                self.assertTrue(any(a.endswith("julia-1.10.12-linux-" + arch + ".tar.gz") for a in downloads[0]["args"]))
                checksum = self.calls("checksum-input")
                self.assertEqual(len(checksum), 1)
                self.assertEqual(checksum[0]["text"].split("  ", 1)[0], expected_hashes[arch])
                self.assertEqual(self.calls("julia")[0]["env"]["SATURN_BOOTSTRAP_ARCH"], arch)
                self.assertTrue((runtime / "julia-1.10.12/bin/julia").is_file())
                self.assertFalse((runtime / ".saturn-julia-bootstrap.lock").exists())
                self.assertFalse(list(runtime.glob(".saturn-julia.*")))
                self.reset_calls()
                self.run_script("slurm/bootstrap.sh", runtime, self.depot, env={"MOCK_ARCH": arch})
                self.assertFalse(self.calls("curl"))
                self.assertEqual(len(self.calls("julia")), 1)

    def test_checksum_failure_prevents_extraction_and_install(self):
        self.run_script("slurm/bootstrap.sh", self.runtime, self.depot,
                        env={"MOCK_BAD_CHECKSUM": "1"}, expected="failure")
        self.assertFalse(self.calls("tar"))
        self.assertFalse(self.calls("julia"))
        self.assertFalse((self.runtime / "julia-1.10.12").exists())
        self.assertFalse((self.runtime / ".saturn-julia-bootstrap.lock").exists())
        self.assertFalse(list(self.runtime.iterdir()))

    def test_unsupported_architecture_stops_before_download(self):
        self.run_script("slurm/bootstrap.sh", self.runtime, self.depot,
                        env={"MOCK_ARCH": "riscv64"}, expected="failure")
        self.assertFalse(self.calls("curl"))
        self.assertFalse(self.calls("julia"))


class ArrayIntegrityTests(SlurmHarness):
    """The common hash-based array contract must survive cluster adapters."""

    def wrappers(self):
        return ("slurm/ensemble.sbatch", "slurm/unity/ensemble.sbatch")

    def test_array_preserves_spaced_config_and_requires_gpu(self):
        manifest, digest = self.manifest()
        for wrapper in self.wrappers():
            with self.subTest(wrapper=wrapper):
                self.reset_calls()
                self.run_script(wrapper, env=self.array_environment(manifest, digest))
                launch = self.assert_gpu_launches()[0]
                run_index = launch["args"].index(str(self.project / "scripts/run.jl"))
                staged = Path(launch["args"][run_index + 1])
                self.assertEqual(staged.read_bytes(), self.spaced_config.read_bytes())
                context = (staged.parent / "job-context.txt").read_text()
                self.assertIn("source_config=" + str(self.spaced_config), context)

    def test_changed_manifest_stops_before_launch(self):
        manifest, digest = self.manifest()
        manifest.chmod(0o600)
        with manifest.open("a") as stream:
            stream.write("tampered\n")
        for wrapper in self.wrappers():
            with self.subTest(wrapper=wrapper):
                self.reset_calls()
                result = self.run_script(wrapper, env=self.array_environment(manifest, digest), expected="failure")
                self.assertIn("manifest has changed", result.stderr)
                self.assertFalse(self.calls("srun"))
                self.assert_gpu_launches(0)

    def test_changed_config_stops_before_launch(self):
        manifest, digest = self.manifest()
        self.spaced_config.write_text("# edited after freezing\n")
        for wrapper in self.wrappers():
            with self.subTest(wrapper=wrapper):
                self.reset_calls()
                result = self.run_script(wrapper, env=self.array_environment(manifest, digest), expected="failure")
                self.assertIn("config has changed", result.stderr)
                self.assertFalse(self.calls("srun"))
                self.assert_gpu_launches(0)

    def test_invalid_array_index_stops_before_launch(self):
        manifest, digest = self.manifest()
        for wrapper in self.wrappers():
            for task in ("0", "3", "not-an-index"):
                with self.subTest(wrapper=wrapper, task=task):
                    self.reset_calls()
                    self.run_script(wrapper, env=self.array_environment(manifest, digest, task), expected="failure")
                    self.assertFalse(self.calls("srun"))
                    self.assert_gpu_launches(0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
