#!/usr/bin/env bash
# A fixed, durable ledger keeps this three-job follow-up separate from other work.
set -euo pipefail
umask 027
source "$(dirname -- "${BASH_SOURCE[0]}")/environment.sh"
[[ $# == 1 && ( $1 == submit || $1 == status ) ]] || unity_die 'usage: followup_01.sh submit | status'
ledger="$SATURN_UNITY_ROOT/followup-01-jobs.tsv"
submission_log="$SATURN_UNITY_ROOT/followup-01-submission.log"

submit_one() {
    local role=$1 mode=$2 config=$3 response job
    shift 3
    if response=$(bash "$SATURN_PROJECT/slurm/unity/submit.sh" "$mode" "$config" "$@" 2>&1); then
        printf '%s\n' "$response"
        printf '%s\n' "$response" >> "$submission_log"
    else
        printf '%s\n' "$response" >&2
        printf '%s\n' "$response" >> "$submission_log"
        unity_die "submission stopped at $role; keep $ledger and $submission_log for recovery; do not resubmit blindly"
    fi
    IFS= read -r job < "$SATURN_UNITY_ROOT/last-job-id"
    [[ $job =~ ^[0-9]+$ ]] || unity_die 'invalid recorded Slurm job ID; inspect the preserved submission log'
    # Also match the response, so a concurrent submission cannot silently replace
    # last-job-id with an unrelated job between the child exit and this read.
    [[ $'\n'"$response"$'\n' == *$'\n'"Submitted Unity $mode job $job"$'\n'* ]] || unity_die 'job-ID record changed unexpectedly; inspect the preserved submission log'
    printf '%s\t%s\n' "$role" "$job" >> "$ledger"
    submitted_job=$job
}

if [[ $1 == submit ]]; then
    [[ -z ${SLURM_JOB_ID:-} ]] || unity_die 'submit from a Unity login shell outside an allocation'
    [[ -z ${SATURN_RESUME:-} ]] || unity_die 'unset SATURN_RESUME; this follow-up starts fresh and creates its own restart checkpoint'
    mkdir -p -- "$SATURN_UNITY_ROOT"
    # Reserve before the first scheduler mutation. Never remove this ledger on
    # failure: even a partial submission must block accidental duplicate jobs.
    ( set -o noclobber; printf 'role\tjob_id\n' > "$ledger" ) 2>/dev/null ||
        unity_die "follow-up already reserved: $ledger; use 'followup_01.sh status' and inspect it before any recovery"
    printf 'Submission ledger: %s\n' "$ledger"
    submit_one restart_check restart-check configs/unity_benchmark.toml
    restart_job=$submitted_job
    submit_one maintenance_forced run configs/unity_followup_01/maintenance_forced.toml "--dependency=afterok:$restart_job"
    submit_one maintenance_unforced run configs/unity_followup_01/maintenance_unforced.toml "--dependency=afterok:$restart_job"
    printf '\nSubmitted follow-up jobs (maintenance waits for restart-check success):\n'
    cat -- "$ledger"
    printf 'Logs: %s/logs\nStatus: bash slurm/unity/followup_01.sh status\n' "$SATURN_RUN_ROOT"
    exit 0
fi

[[ -f $ledger ]] || unity_die "no follow-up ledger found: $ledger"
IFS= read -r header < "$ledger"
[[ $header == $'role\tjob_id' ]] || unity_die "invalid ledger header: $ledger"
job_list= seen='|'
while IFS=$'\t' read -r role job extra; do
    [[ $role != role ]] || continue
    case $role in restart_check|maintenance_forced|maintenance_unforced) ;; *) unity_die "invalid role in $ledger" ;; esac
    [[ $job =~ ^[0-9]+$ && -z $extra ]] || unity_die "invalid job row in $ledger"
    [[ $seen != *"|$role|"* ]] || unity_die "duplicate role in $ledger"
    seen+="$role|"
    job_list="${job_list:+$job_list,}$job"
    printf '%s: %s\n' "$role" "$job"
done < "$ledger"
printf 'Ledger: %s\nLogs: %s/logs\n' "$ledger" "$SATURN_RUN_ROOT"
[[ -n $job_list ]] || { printf 'No confirmed job IDs were recorded; inspect %s before recovery.\n' "$submission_log"; exit 0; }
[[ $seen == '|restart_check|maintenance_forced|maintenance_unforced|' ]] ||
    printf 'Partial submission: only the recorded jobs below are being queried.\n'
# Finished jobs may no longer be visible to squeue; still query accounting.
queue_status=0
squeue -j "$job_list" -o '%.18i %.30j %.12T %.10M %R' || queue_status=$?
sacct -j "$job_list" -X --format=JobID,JobName,State,ExitCode,Elapsed,Reason
(( queue_status == 0 )) || printf 'squeue did not find all jobs; use the accounting rows above for completed jobs.\n'
