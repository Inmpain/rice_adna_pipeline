#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SUT="${SCRIPT_DIR}/submit_oryza_competitive_mapping.sh"
TEST_ROOT="$(mktemp -d /tmp/oryza_mapping_test.XXXXXX)"
trap 'rm -rf "$TEST_ROOT"' EXIT

mkdir -p \
    "$TEST_ROOT/reads" \
    "$TEST_ROOT/wgs" \
    "$TEST_ROOT/rice" \
    "$TEST_ROOT/out"

for shard in $(seq 1 129); do
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

sbatch() {
    case " $* " in
        *" --test-only "*) return 0 ;;
        *) printf '9001\n' ;;
    esac
}
sinfo() { printf 'comp*\n'; }
bowtie2() { return 0; }
samtools() { return 0; }
export -f sbatch sinfo bowtie2 samtools

COMMON_ENV=(
    READ_DIR="$TEST_ROOT/reads"
    WGS_DB_DIR="$TEST_ROOT/wgs"
    RICE_DB_DIR="$TEST_ROOT/rice"
    OUT_DIR="$TEST_ROOT/out"
    SLURM_PARTITION=comp
    SLURM_ACCOUNT=
)

env "${COMMON_ENV[@]}" bash "$SUT" check > "$TEST_ROOT/check.log"
grep -Fq '[check] samples=1 databases=131' "$TEST_ROOT/check.log"
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

set +e
env "${COMMON_ENV[@]}" bash "$SUT" submit > "$TEST_ROOT/locked.log" 2>&1
locked_status="$?"
set -e
[[ "$locked_status" -eq 2 ]]
grep -Fq 'full submission is locked' "$TEST_ROOT/locked.log"
printf 'PASS full_submit_lock\n'

FULL_OUT="$TEST_ROOT/full_out"
env \
    READ_DIR="$TEST_ROOT/reads" \
    WGS_DB_DIR="$TEST_ROOT/wgs" \
    RICE_DB_DIR="$TEST_ROOT/rice" \
    OUT_DIR="$FULL_OUT" \
    SLURM_PARTITION=comp \
    SLURM_ACCOUNT= \
    CONFIRM_FULL_SUBMIT=YES \
    bash "$SUT" submit > "$TEST_ROOT/full_submit.log"
submission_file="$(find "$FULL_OUT/submissions" -name 'submitted_jobs.*.tsv' -print -quit)"
[[ -n "$submission_file" ]]
[[ "$(wc -l < "$submission_file" | awk '{print $1}')" -eq 133 ]]
grep -Fq $'merge_sort\tS1\tALL\t9001' "$submission_file"
printf 'PASS confirmed_full_submit_plan_131_maps_plus_merge\n'
