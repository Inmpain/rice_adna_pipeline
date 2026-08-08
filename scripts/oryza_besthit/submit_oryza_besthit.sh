#!/usr/bin/env bash
#
# Submit the Oryza best-hit / aDNA-damage-aware filtering stage on SLURM.
# One SLURM job per sample; each job runs oryza_besthit_damage_filter.py once
# against that sample's query-name sorted competitive-mapping BAM.
#
# check mode validates every input path, python3 + pysam, and the SLURM
# partition/account, without submitting anything.
# smoke mode submits exactly one job for one sample with --limit-reads 1000,
# writing to an isolated smoke_test output dir (never touches production).
# submit mode submits one job per given sample (skips samples that already
# have a <sample>.finished marker in OUT_DIR, so re-running is safe).
# merge mode concatenates every <sample>.summary.tsv in OUT_DIR into
# besthit_summary.tsv. Run it only after all submit jobs have finished --
# it deliberately does NOT run inside a job, so parallel per-sample jobs never
# race on one shared summary file.
#
# The script is also its own SLURM worker (internal `run` mode).

set -Eeuo pipefail

report_error() {
    local exit_code="$?"
    echo "ERROR: line ${BASH_LINENO[0]} failed with exit ${exit_code}: ${BASH_COMMAND}" >&2
    exit "$exit_code"
}
trap report_error ERR

# -----------------------------------------------------------------------------
# Paths. Environment variables with the same names override these defaults.
# -----------------------------------------------------------------------------

BAM_DIR="${BAM_DIR:-/home/scratch/yinmt202607/gene/results/oryza_competitive_mapping/by_sample}"
BAM_SUFFIX="${BAM_SUFFIX:-.competitive.name_sorted.bam}"

FASTQ_DIR="${FASTQ_DIR:-/home/scratch/yinmt202607/gene/results/oryza_candidates_combined}"
FASTQ_SUFFIX="${FASTQ_SUFFIX:-.oryza_candidates.combined.fastq.gz}"

ACC2TAXID="${ACC2TAXID:-/home/scratch/yinmt202607/db/asian_rice_panel_index/all_wgs_asian_irgsp.acc2taxid}"
NODES="${NODES:-/home/database/ref20250728/taxonomy_CPH/ncbi/20250530/nodes.dmp}"
NAMES="${NAMES:-/home/database/ref20250728/taxonomy_CPH/ncbi/20250530/names.dmp}"

# v2 (2026-08-08): Oryza scope defaults to the WHOLE genus now, not a
# hardcoded 3-species list -- see oryza_besthit_damage_filter.py's module
# docstring and docs/ORYZA_BESTHIT_HANDOFF.md section 5.1b. Leave
# ORYZA_TAXIDS empty (the default) for genus-wide auto-resolution driven by
# ORYZA_GENUS_TAXID; set ORYZA_TAXIDS explicitly (e.g. "4529 4530 4536") to
# reproduce v1's narrower rufipogon/sativa/nivara-only behavior instead.
ORYZA_TAXIDS="${ORYZA_TAXIDS:-}"
ORYZA_GENUS_TAXID="${ORYZA_GENUS_TAXID:-4527}"
DAMAGE_WINDOW="${DAMAGE_WINDOW:-5}"
TOP_N="${TOP_N:-10}"
# Optional pre-gate, off by default (both empty) -- see --min-best-similarity/
# --max-best-raw-nm in the python script's --help.
MIN_BEST_SIMILARITY="${MIN_BEST_SIMILARITY:-}"
MAX_BEST_RAW_NM="${MAX_BEST_RAW_NM:-}"

OUT_DIR="${OUT_DIR:-/home/scratch/yinmt202607/gene/results/oryza_competitive_mapping/besthit}"
LOG_DIR="${OUT_DIR}/logs"
SUBMIT_DIR="${OUT_DIR}/submissions"

SLURM_PARTITION="${SLURM_PARTITION:-comp}"
SLURM_ACCOUNT="${SLURM_ACCOUNT:-}"

JOB_CPUS="${JOB_CPUS:-4}"
JOB_MEM_MB="${JOB_MEM_MB:-16000}"
JOB_TIME="${JOB_TIME:-04:00:00}"

TEST_READS="${TEST_READS:-1000}"
TEST_OUT_DIR="${TEST_OUT_DIR:-${OUT_DIR}/smoke_test}"
TEST_TIME="${TEST_TIME:-00:30:00}"

SCRIPT_PATH="$(readlink -f "${BASH_SOURCE[0]}")"
SCRIPT_DIR="$(dirname "$SCRIPT_PATH")"
PY_SCRIPT="${PY_SCRIPT:-${SCRIPT_DIR}/oryza_besthit_damage_filter.py}"

usage() {
    cat <<'EOF'
Usage:
  submit_oryza_besthit.sh check
      Validate inputs, python3+pysam, and the SLURM partition/account.
      Submits nothing.

  submit_oryza_besthit.sh smoke [SAMPLE]
      Submit one job for one sample (first discovered if omitted) with
      --limit-reads 1000. Output goes to OUT_DIR/smoke_test, isolated from
      production. Does not write a .finished marker (see the python script).

  submit_oryza_besthit.sh submit SAMPLE [SAMPLE ...]
      Submit one full job per sample. Skips a sample if
      OUT_DIR/<sample>.finished already exists.

  submit_oryza_besthit.sh submit all
      Same, but auto-discovers every sample with a finished competitive-
      mapping BAM under BAM_DIR (same discovery `check` uses). Safe to
      re-run as more samples finish mapping -- already-finished besthit
      samples are skipped.

  submit_oryza_besthit.sh local SAMPLE [SAMPLE ...] | local all
      Same as submit, but runs sequentially in the foreground -- no sbatch,
      no queue wait. Good for quick iteration/debugging. Aborts on the
      first failing sample. Log still goes to
      OUT_DIR/logs/<sample>.besthit.local.log (as well as the terminal).

  submit_oryza_besthit.sh merge
      Concatenate OUT_DIR/<sample>.summary.tsv (all samples that have run,
      finished or not) into OUT_DIR/besthit_summary.tsv. Run after jobs land.
EOF
}

report_error_reset() { trap report_error ERR; }

timestamp() {
    date --iso-8601=seconds 2>/dev/null || date '+%Y-%m-%dT%H:%M:%S%z'
}

load_python_module() {
    if command -v module >/dev/null 2>&1; then
        module load python/ >/dev/null 2>&1 || true
    fi
    command -v python3 >/dev/null 2>&1 || {
        echo "ERROR: python3 not found on PATH" >&2
        exit 1
    }
}

check_pysam() {
    python3 -c "import pysam" >/dev/null 2>&1 || {
        echo "ERROR: 'python3 -c \"import pysam\"' failed. Install pysam in" >&2
        echo "  the environment this script's python3 resolves to (e.g." >&2
        echo "  'pip install --user pysam' or a dedicated conda env) before" >&2
        echo "  submitting." >&2
        exit 1
    }
}

sbatch_common_args() {
    SBATCH_COMMON=(--partition="$SLURM_PARTITION")
    if [[ -n "$SLURM_ACCOUNT" ]]; then
        SBATCH_COMMON+=(--account="$SLURM_ACCOUNT")
    fi
}

bam_path() {
    printf '%s/%s%s\n' "$BAM_DIR" "$1" "$BAM_SUFFIX"
}

fastq_path() {
    printf '%s/%s%s\n' "$FASTQ_DIR" "$1" "$FASTQ_SUFFIX"
}

discover_samples() {
    [[ -d "$BAM_DIR" ]] || {
        echo "ERROR: BAM_DIR not found: $BAM_DIR" >&2
        return 1
    }
    shopt -s nullglob
    local finished=("$BAM_DIR"/*"${BAM_SUFFIX}.finished")
    shopt -u nullglob
    SAMPLES=()
    local f base
    for f in "${finished[@]}"; do
        base="$(basename "$f")"
        SAMPLES+=("${base%${BAM_SUFFIX}.finished}")
    done
}

validate_slurm_request() {
    command -v sbatch >/dev/null 2>&1 || {
        echo "ERROR: sbatch not found on PATH" >&2
        return 1
    }
    sbatch_common_args
    if command -v sinfo >/dev/null 2>&1; then
        local available
        available="$(sinfo -h -o '%P' | tr -d '*' | sort -u)"
        if ! grep -Fxq "$SLURM_PARTITION" <<< "$available"; then
            echo "ERROR: SLURM partition does not exist: $SLURM_PARTITION" >&2
            echo "Available partitions:" >&2
            sed 's/^/  /' <<< "$available" >&2
            return 1
        fi
    fi
    sbatch --test-only "${SBATCH_COMMON[@]}" \
        --job-name=orybh.preflight --cpus-per-task="$JOB_CPUS" \
        --mem="${JOB_MEM_MB}M" --time="$JOB_TIME" --wrap=true >/dev/null
}

run_check() {
    echo "[check] BAM_DIR=$BAM_DIR"
    echo "[check] FASTQ_DIR=$FASTQ_DIR"
    echo "[check] ACC2TAXID=$ACC2TAXID"
    echo "[check] NODES=$NODES"
    echo "[check] NAMES=$NAMES"
    if [[ -n "$ORYZA_TAXIDS" ]]; then
        echo "[check] Oryza scope: manual override ORYZA_TAXIDS=$ORYZA_TAXIDS"
    else
        echo "[check] Oryza scope: whole genus, ORYZA_GENUS_TAXID=$ORYZA_GENUS_TAXID"
    fi
    echo "[check] DAMAGE_WINDOW=$DAMAGE_WINDOW TOP_N=$TOP_N"
    if [[ -n "$MIN_BEST_SIMILARITY" || -n "$MAX_BEST_RAW_NM" ]]; then
        echo "[check] quality pre-gate ON: MIN_BEST_SIMILARITY=$MIN_BEST_SIMILARITY MAX_BEST_RAW_NM=$MAX_BEST_RAW_NM"
    fi
    echo "[check] partition=$SLURM_PARTITION account=${SLURM_ACCOUNT:-<default>}"

    [[ -f "$ACC2TAXID" ]] || { echo "ERROR: ACC2TAXID missing: $ACC2TAXID" >&2; exit 1; }
    [[ -f "$NODES" ]] || { echo "ERROR: NODES missing: $NODES" >&2; exit 1; }
    [[ -f "$NAMES" ]] || { echo "ERROR: NAMES missing: $NAMES" >&2; exit 1; }
    [[ -f "$PY_SCRIPT" ]] || { echo "ERROR: PY_SCRIPT missing: $PY_SCRIPT" >&2; exit 1; }

    discover_samples
    [[ "${#SAMPLES[@]}" -gt 0 ]] || {
        echo "ERROR: no *${BAM_SUFFIX}.finished markers under $BAM_DIR" >&2
        exit 1
    }
    echo "[check] samples with a finished competitive-mapping BAM: ${#SAMPLES[@]}"
    printf '  %s\n' "${SAMPLES[@]}"

    local s
    for s in "${SAMPLES[@]}"; do
        [[ -f "$(fastq_path "$s")" ]] || {
            echo "ERROR: candidate FASTQ missing for $s: $(fastq_path "$s")" >&2
            exit 1
        }
    done

    load_python_module
    check_pysam
    validate_slurm_request
    echo "[check] PASS"
}

run_worker() {
    # Internal: invoked as the sbatch job body. Args: sample outdir [--limit-reads N]
    [[ "$#" -ge 2 ]] || { usage >&2; exit 2; }
    local sample="$1" outdir="$2"; shift 2
    local bam fastq
    bam="$(bam_path "$sample")"
    fastq="$(fastq_path "$sample")"

    load_python_module
    mkdir -p "$outdir"

    echo "[run] sample=$sample started=$(timestamp)"
    echo "[run] bam=$bam"
    echo "[run] fastq=$fastq"
    echo "[run] outdir=$outdir"

    # v2: pass an explicit whitelist only if ORYZA_TAXIDS is set; otherwise
    # let the python script auto-resolve the whole genus via
    # --oryza-genus-taxid (its default, 4527, is already genus Oryza -- we
    # still pass ORYZA_GENUS_TAXID explicitly so a non-default env var value
    # actually takes effect).
    local oryza_args=()
    if [[ -n "$ORYZA_TAXIDS" ]]; then
        # shellcheck disable=SC2206
        oryza_args=(--oryza-taxids $ORYZA_TAXIDS)
    else
        oryza_args=(--oryza-genus-taxid "$ORYZA_GENUS_TAXID")
    fi
    local gate_args=()
    [[ -n "$MIN_BEST_SIMILARITY" ]] && gate_args+=(--min-best-similarity "$MIN_BEST_SIMILARITY")
    [[ -n "$MAX_BEST_RAW_NM" ]] && gate_args+=(--max-best-raw-nm "$MAX_BEST_RAW_NM")

    python3 "$PY_SCRIPT" \
        --sample "$sample" \
        --bam "$bam" \
        --fastq "$fastq" \
        --acc2taxid "$ACC2TAXID" \
        --nodes "$NODES" \
        --names "$NAMES" \
        "${oryza_args[@]}" \
        --damage-window "$DAMAGE_WINDOW" \
        --top-n "$TOP_N" \
        "${gate_args[@]}" \
        --outdir "$outdir" \
        --threads "$JOB_CPUS" \
        "$@"

    echo "[run] sample=$sample finished=$(timestamp)"
}

run_smoke() {
    load_python_module
    check_pysam
    discover_samples
    local sample="${1:-}"
    if [[ -z "$sample" ]]; then
        [[ "${#SAMPLES[@]}" -gt 0 ]] || {
            echo "ERROR: no finished samples to pick a default from" >&2
            exit 1
        }
        sample="${SAMPLES[0]}"
    fi
    [[ -f "$(bam_path "$sample")" ]] || {
        echo "ERROR: BAM missing for $sample: $(bam_path "$sample")" >&2
        exit 1
    }

    mkdir -p "$TEST_OUT_DIR" "$LOG_DIR" "$SUBMIT_DIR"
    validate_slurm_request
    sbatch_common_args

    local log="${LOG_DIR}/${sample}.smoke.%j.log"
    local job_id
    job_id="$(sbatch \
        --parsable \
        "${SBATCH_COMMON[@]}" \
        --export="ALL,PY_SCRIPT=${PY_SCRIPT}" \
        --job-name="orybh.smoke.${sample}" \
        --cpus-per-task="$JOB_CPUS" \
        --mem="${JOB_MEM_MB}M" \
        --time="$TEST_TIME" \
        --output="$log" \
        --error="$log" \
        "$SCRIPT_PATH" run "$sample" "$TEST_OUT_DIR" --limit-reads "$TEST_READS")"

    echo "[smoke] submitted job_id=$job_id sample=$sample reads=$TEST_READS"
    echo "[smoke] watch: squeue -j $job_id"
    echo "[smoke] log:   tail -f ${log//%j/$job_id}"
    echo "[smoke] output: $TEST_OUT_DIR (no .finished marker written; expected)"
}

run_submit() {
    [[ "$#" -ge 1 ]] || { usage >&2; exit 2; }
    load_python_module
    check_pysam
    mkdir -p "$OUT_DIR" "$LOG_DIR" "$SUBMIT_DIR"
    validate_slurm_request
    sbatch_common_args

    local manifest="${SUBMIT_DIR}/submit.$(date +%Y%m%dT%H%M%S).tsv"
    printf 'job_id\tsample\n' > "$manifest"

    local sample bam fastq log job_id
    for sample in "$@"; do
        bam="$(bam_path "$sample")"
        fastq="$(fastq_path "$sample")"
        [[ -f "$bam" ]] || { echo "ERROR: BAM missing for $sample: $bam" >&2; exit 1; }
        [[ -f "$fastq" ]] || { echo "ERROR: FASTQ missing for $sample: $fastq" >&2; exit 1; }

        if [[ -f "${OUT_DIR}/${sample}.finished" ]]; then
            echo "[submit] $sample already finished, skipping"
            continue
        fi

        log="${LOG_DIR}/${sample}.besthit.%j.log"
        job_id="$(sbatch \
            --parsable \
            "${SBATCH_COMMON[@]}" \
            --export="ALL,PY_SCRIPT=${PY_SCRIPT}" \
            --job-name="orybh.${sample}" \
            --cpus-per-task="$JOB_CPUS" \
            --mem="${JOB_MEM_MB}M" \
            --time="$JOB_TIME" \
            --output="$log" \
            --error="$log" \
            "$SCRIPT_PATH" run "$sample" "$OUT_DIR")"
        printf '%s\t%s\n' "$job_id" "$sample" >> "$manifest"
        echo "[submit] sample=$sample job_id=$job_id"
    done
    echo "[submit] manifest: $manifest"
}

run_local() {
    # Sequential, foreground, no sbatch -- same worker + .finished skip logic
    # as run_submit, just without a SLURM job in between. Useful for quick
    # iteration/debugging on a login node or an already-allocated interactive
    # session. Aborts on the first failing sample (same fail-fast convention
    # as the rest of this script) rather than silently skipping ahead.
    [[ "$#" -ge 1 ]] || { usage >&2; exit 2; }
    load_python_module
    check_pysam
    mkdir -p "$OUT_DIR" "$LOG_DIR"

    local sample bam fastq log
    for sample in "$@"; do
        bam="$(bam_path "$sample")"
        fastq="$(fastq_path "$sample")"
        [[ -f "$bam" ]] || { echo "ERROR: BAM missing for $sample: $bam" >&2; exit 1; }
        [[ -f "$fastq" ]] || { echo "ERROR: FASTQ missing for $sample: $fastq" >&2; exit 1; }

        if [[ -f "${OUT_DIR}/${sample}.finished" ]]; then
            echo "[local] $sample already finished, skipping"
            continue
        fi

        log="${LOG_DIR}/${sample}.besthit.local.log"
        echo "[local] running sample=$sample (foreground, no sbatch) -> $log"
        run_worker "$sample" "$OUT_DIR" 2>&1 | tee "$log"
    done
}

run_merge() {
    shopt -s nullglob
    local files=("$OUT_DIR"/*.summary.tsv)
    shopt -u nullglob
    [[ "${#files[@]}" -gt 0 ]] || {
        echo "ERROR: no *.summary.tsv under $OUT_DIR yet" >&2
        exit 1
    }

    local out="${OUT_DIR}/besthit_summary.tsv"
    local tmp="${out}.tmp.$$"
    rm -f "$tmp"
    trap 'rm -f "$tmp"' EXIT

    head -n1 "${files[0]}" > "$tmp"
    local f
    for f in "${files[@]}"; do
        tail -n +2 "$f" >> "$tmp"
    done
    mv "$tmp" "$out"
    trap - EXIT
    report_error_reset
    echo "[merge] ${#files[@]} samples -> $out"
}

case "${1:-}" in
    check)
        run_check
        ;;
    smoke)
        run_smoke "${2:-}"
        ;;
    submit)
        shift
        if [[ "${1:-}" == "all" ]]; then
            discover_samples
            [[ "${#SAMPLES[@]}" -gt 0 ]] || {
                echo "ERROR: no *${BAM_SUFFIX}.finished markers under $BAM_DIR" >&2
                exit 1
            }
            echo "[submit] all: ${#SAMPLES[@]} samples with a finished mapping BAM"
            run_submit "${SAMPLES[@]}"
        else
            run_submit "$@"
        fi
        ;;
    local)
        shift
        if [[ "${1:-}" == "all" ]]; then
            discover_samples
            [[ "${#SAMPLES[@]}" -gt 0 ]] || {
                echo "ERROR: no *${BAM_SUFFIX}.finished markers under $BAM_DIR" >&2
                exit 1
            }
            echo "[local] all: ${#SAMPLES[@]} samples with a finished mapping BAM"
            run_local "${SAMPLES[@]}"
        else
            run_local "$@"
        fi
        ;;
    merge)
        run_merge
        ;;
    run)
        shift
        run_worker "$@"
        ;;
    *)
        usage
        exit 2
        ;;
esac
