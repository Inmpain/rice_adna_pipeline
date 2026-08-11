#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SUT="${SCRIPT_DIR}/submit_oryza_competitive_mapping.sh"
TEST_ROOT="$(mktemp -d /tmp/oryza_mapping_test.XXXXXX)"
trap 'rm -rf "$TEST_ROOT"' EXIT
export WGS_FIRST_SHARD=1
export WGS_LAST_SHARD=2

mkdir -p \
    "$TEST_ROOT/reads" \
    "$TEST_ROOT/wgs" \
    "$TEST_ROOT/rice" \
    "$TEST_ROOT/out"

for shard in $(seq "$WGS_FIRST_SHARD" "$WGS_LAST_SHARD"); do
    prefix="$TEST_ROOT/wgs/wgs_eukaryota.${shard}.fas.gz"
    for suffix in .1.bt2l .2.bt2l .3.bt2l .4.bt2l .rev.1.bt2l .rev.2.bt2l; do
        truncate -s 1024 "${prefix}${suffix}"
    done
done

for prefix in \
    "$TEST_ROOT/rice/asian_rice_panel.fa" \
    "$TEST_ROOT/rice/irgsp_bt2idx"
do
    for suffix in .1.bt2l .2.bt2l .3.bt2l .4.bt2l .rev.1.bt2l .rev.2.bt2l; do
        truncate -s 1024 "${prefix}${suffix}"
    done
done

printf '@r1\nACGT\n+\nIIII\n@r2\nTGCA\n+\nIIII\n' \
    | gzip -c \
    > "$TEST_ROOT/reads/S1.oryza_candidates.combined.fastq.gz"
printf '@r3\nAAAA\n+\nIIII\n@r4\nCCCC\n+\nIIII\n' \
    | gzip -c \
    > "$TEST_ROOT/reads/S2.oryza_candidates.combined.fastq.gz"

sbatch() {
    if [[ -n "${FAKE_SBATCH_LOG:-}" ]]; then
        printf '%s\n' "$*" >> "$FAKE_SBATCH_LOG"
    fi
    case " $* " in
        *" --test-only "*) return 0 ;;
        *) printf '9001\n' ;;
    esac
}
sinfo() { printf 'comp*\n'; }
squeue() {
    if [[ "${FAKE_ACTIVE:-0}" == "1" ]]; then
        printf '777|orymap.S1.wgs_eukaryota.1|RUNNING\n'
    fi
}
bowtie2() { return 0; }
samtools() { return 0; }
export -f sbatch sinfo squeue bowtie2 samtools

COMMON_ENV=(
    READ_DIR="$TEST_ROOT/reads"
    WGS_DB_DIR="$TEST_ROOT/wgs"
    RICE_DB_DIR="$TEST_ROOT/rice"
    OUT_DIR="$TEST_ROOT/out"
    SLURM_PARTITION=comp
    SLURM_ACCOUNT=
    FAKE_SBATCH_LOG="$TEST_ROOT/fake_sbatch.log"
)

env "${COMMON_ENV[@]}" bash "$SUT" check > "$TEST_ROOT/check.log"
grep -Fq '[check] samples=2 databases=4' "$TEST_ROOT/check.log"
grep -Fq '[check] PASS:' "$TEST_ROOT/check.log"
printf 'PASS mock_check\n'

if env "${COMMON_ENV[@]}" SLURM_PARTITION=wrong \
    bash "$SUT" check > "$TEST_ROOT/wrong.log" 2>&1
then
    echo 'ERROR: wrong partition was not rejected' >&2
    exit 1
fi
grep -Fq 'SLURM partition does not exist: wrong' "$TEST_ROOT/wrong.log"
printf 'PASS invalid_partition_guard\n'

env "${COMMON_ENV[@]}" TEST_READS=1 \
    bash "$SUT" test > "$TEST_ROOT/smoke.log"
grep -Fq '[test] submitted exactly one mapping job' "$TEST_ROOT/smoke.log"
grep -Fq 'job_id=9001 database=irgsp' "$TEST_ROOT/smoke.log"
[[ -s "$TEST_ROOT/out/smoke_test/input/S1.smoke_test.first_1_reads.fastq.gz" ]]
printf 'PASS one_job_smoke_test\n'

[[ "$(env "${COMMON_ENV[@]}" bash "$SUT" list)" == $'S1\nS2' ]]
printf 'PASS list_samples\n'

set +e
env "${COMMON_ENV[@]}" bash "$SUT" submit > "$TEST_ROOT/no_sample.log" 2>&1
no_sample_status="$?"
set -e
[[ "$no_sample_status" -eq 2 ]]
grep -Fq 'submit LV6000619499' "$TEST_ROOT/no_sample.log"
printf 'PASS all_sample_submit_removed\n'

set +e
env "${COMMON_ENV[@]}" bash "$SUT" submit UNKNOWN \
    > "$TEST_ROOT/unknown_sample.log" 2>&1
unknown_status="$?"
set -e
[[ "$unknown_status" -eq 1 ]]
grep -Fq 'sample FASTQ not found' "$TEST_ROOT/unknown_sample.log"
printf 'PASS unknown_sample_rejected\n'

set +e
env "${COMMON_ENV[@]}" FAKE_ACTIVE=1 bash "$SUT" submit S1 \
    > "$TEST_ROOT/active_sample.log" 2>&1
active_status="$?"
set -e
[[ "$active_status" -ne 0 ]]
grep -Fq 'already has queued/running workflow jobs' "$TEST_ROOT/active_sample.log"
printf 'PASS duplicate_active_sample_guard\n'

FULL_OUT="$TEST_ROOT/full_out"
env \
    READ_DIR="$TEST_ROOT/reads" \
    WGS_DB_DIR="$TEST_ROOT/wgs" \
    RICE_DB_DIR="$TEST_ROOT/rice" \
    OUT_DIR="$FULL_OUT" \
    SLURM_PARTITION=comp \
    SLURM_ACCOUNT= \
    bash "$SUT" submit S1 > "$TEST_ROOT/sample_submit.log"
submission_file="$(find "$FULL_OUT/submissions" -name 'submitted_jobs.*.tsv' -print -quit)"
[[ -n "$submission_file" ]]
[[ "$(wc -l < "$submission_file" | awk '{print $1}')" -eq 6 ]]
grep -Fq $'merge_sort\tS1\tALL\t9001' "$submission_file"
grep -Fq '[submit] sample=S1 databases=4' "$TEST_ROOT/sample_submit.log"
grep -Fq 'this invocation will submit at most 5 jobs' "$TEST_ROOT/sample_submit.log"
printf 'PASS one_sample_maps_plus_merge\n'

set +e
env "${COMMON_ENV[@]}" FAKE_ACTIVE=1 bash "$SUT" series --all \
    > "$TEST_ROOT/active_series.log" 2>&1
active_series_status="$?"
set -e
[[ "$active_series_status" -ne 0 ]]
grep -Fq 'production workflow jobs are still queued/running' "$TEST_ROOT/active_series.log"
printf 'PASS active_workflow_blocks_new_series\n'

SERIES_OUT="$TEST_ROOT/series_out"
: > "$TEST_ROOT/series_sbatch.log"
env \
    READ_DIR="$TEST_ROOT/reads" \
    WGS_DB_DIR="$TEST_ROOT/wgs" \
    RICE_DB_DIR="$TEST_ROOT/rice" \
    OUT_DIR="$SERIES_OUT" \
    SLURM_PARTITION=comp \
    SLURM_ACCOUNT= \
    FAKE_SBATCH_LOG="$TEST_ROOT/series_sbatch.log" \
    bash "$SUT" series --all > "$TEST_ROOT/series_start.log"

grep -Fq '[series] id=' "$TEST_ROOT/series_start.log"
grep -Fq 'only the first incomplete sample will be submitted now' "$TEST_ROOT/series_start.log"
grep -Fq 'next sample=S2; later samples are not submitted yet' "$TEST_ROOT/series_start.log"
grep -Fq -- '--job-name=orynext.S2' "$TEST_ROOT/series_sbatch.log"
grep -Fq -- 'series-next' "$TEST_ROOT/series_sbatch.log"
[[ "$(find "$SERIES_OUT/submissions" -name 'submitted_jobs.S1.*.tsv' | wc -l | awk '{print $1}')" -eq 1 ]]
[[ "$(find "$SERIES_OUT/submissions" -name 'submitted_jobs.S2.*.tsv' | wc -l | awk '{print $1}')" -eq 0 ]]

series_state="$(find "$SERIES_OUT/series" -name 'series.*.state.tsv' -print -quit)"
[[ -n "$series_state" ]]
grep -Fq $'continuation_submitted\tS1\t9001\t9001\tS2\t1' "$series_state"
series_id="$(basename "$series_state")"
series_id="${series_id#series.}"
series_id="${series_id%.state.tsv}"
printf 'PASS series_starts_only_first_sample\n'

env \
    READ_DIR="$TEST_ROOT/reads" \
    WGS_DB_DIR="$TEST_ROOT/wgs" \
    RICE_DB_DIR="$TEST_ROOT/rice" \
    OUT_DIR="$SERIES_OUT" \
    SLURM_PARTITION=comp \
    SLURM_ACCOUNT= \
    FAKE_SBATCH_LOG="$TEST_ROOT/series_sbatch.log" \
    bash "$SUT" series-next "$series_id" S2 > "$TEST_ROOT/series_next.log"

grep -Fq 'final sample submitted; no continuation job is needed' "$TEST_ROOT/series_next.log"
[[ "$(find "$SERIES_OUT/submissions" -name 'submitted_jobs.S2.*.tsv' | wc -l | awk '{print $1}')" -eq 1 ]]
grep -Fq $'final_sample_waiting_for_merge\tS2\t9001' "$series_state"
printf 'PASS continuation_submits_next_sample_only\n'
