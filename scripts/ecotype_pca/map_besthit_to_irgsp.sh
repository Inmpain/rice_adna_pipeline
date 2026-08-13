#!/usr/bin/env bash
#
# Maps besthit-filtered ancient Oryza reads to irgsp.fa -- the shared
# reference coordinate system all three ecotype-pca panels (29M_3k,
# 6.7M_720, Civan) are confirmed to use (see ECOTYPE_PCA_PANEL.md 3.1).
# One BAM per sample feeds pseudo-haploid genotype calling against any
# of the three panels -- mapping is only done once per sample, not once
# per panel.
#
# Input: <sample>.besthit_oryza.fastq.gz (codex/oryza-competitive-mapping
#   branch's ORYZA_BESTHIT_HANDOFF.md section 5.3 output). This FASTQ has
#   NOT been deduplicated yet (explicitly deferred "further downstream"
#   per that doc's section 5.2) -- deduplication happens here, since this
#   is that downstream step.
#
# 2026-08-12 hardening (see docs/ECOTYPE_PCA_EXECUTION_PLAN.md P0-2): the
# mapping/filtering/dedup chain is aligned with the project's main pipeline
# (scripts/server_originals/mapping.sh) instead of a from-literature
# reconstruction that had drifted from it:
#   - bwa aln -l 1024 -n 0.01 -o 2 is UNCHANGED -- this already matched
#     mapping.sh exactly, confirmed by direct diff.
#   - samtools view -bh -F 0x904 is applied right after samse, same
#     position as mapping.sh -- drops unmapped/secondary/supplementary
#     alignments before sorting.
#   - dedup goes through the same collate -> fixmate -m -> sort -> markdup
#     chain as mapping.sh, instead of running markdup directly on a
#     coordinate-sorted BAM.
#   - DELIBERATE remaining difference from mapping.sh: markdup here does
#     NOT use -r (does not remove duplicates), only flags them --
#     pseudo_haploid_call.py already filters flagged duplicates at pileup
#     time (aln.is_duplicate check), and the Phase 0 IRGSP coverage census
#     wants duplicate-rate visible in the BAM as a QC signal.
#   - the per-sample summary line reports mapped-primary reads via
#     samtools view -c -F 4, plus separate MAPQ>=30/MAPQ>=20 counts.
#
# 2026-08-12 SECOND revision: restructured to match this repo's
# established `submit_oryza_besthit.sh` pattern (codex/oryza-competitive-
# mapping branch) -- one SLURM job PER SAMPLE instead of one sequential
# for-loop processing all samples in a single shell/job. Running 16
# samples' bwa aln + dedup chain serially in one process (the previous
# version's behavior) needlessly serializes work that is fully
# sample-independent. Subcommands, mirroring submit_oryza_besthit.sh:
#   check              validate paths/tools/SLURM partition, submit nothing
#   smoke [SAMPLE]     submit exactly one job for one sample
#   submit SAMPLE...   one sbatch job per given sample
#   submit all         auto-discover every sample with a besthit FASTQ
#                       under BESTHIT_DIR and submit one job each
#   local SAMPLE...|all  same per-sample worker, sequential, no sbatch --
#                       for quick debugging on an interactive/login shell
#   run SAMPLE OUTDIR  internal: the actual per-sample worker (sbatch job
#                       body); do not call this directly
# Both submit and local skip a sample if OUT_DIR/<sample>.finished already
# exists, so re-running after a partial batch is safe.
#
# Tools: bwa is expected already on PATH (base conda env on this cluster,
# confirmed 2026-08-12 -- no module-load convention for bwa exists
# anywhere in this repo). samtools needs `module load samtools`, which
# load_tools() below does automatically, same convention as
# scripts/oryza_besthit/submit_oryza_competitive_mapping.sh's
# load_samtools_module() (codex/oryza-competitive-mapping branch).
#
# 2026-08-13 THIRD revision: generalized to support the ORSC-narrowed
# target read set (see scripts/oryza_besthit/split_besthit_taxonomic_tiers.py)
# without touching the already-finished whole-genus BAMs. INPUT_SUFFIX and
# READSET_LABEL default to the exact previous hardcoded values
# (.besthit_oryza.fastq.gz / besthit_oryza), so default invocation is
# byte-for-byte unchanged -- point BESTHIT_DIR/INPUT_SUFFIX/READSET_LABEL/
# OUT_DIR at the taxonomic_tiers output + a NEW out dir to run the ORSC
# readset instead (see docs/ECOTYPE_PCA_PHASE0_COMMANDS.md). Two readsets
# must use different OUT_DIR -- the .finished marker is not namespaced by
# READSET_LABEL, so sharing OUT_DIR between two readsets for the same
# sample would make the second readset's run silently skip.

set -Eeuo pipefail

report_error() {
    local exit_code="$?"
    echo "ERROR: line ${BASH_LINENO[0]} failed with exit ${exit_code}: ${BASH_COMMAND}" >&2
    exit "$exit_code"
}
trap report_error ERR

# -----------------------------------------------------------------------------
# Paths / config. Environment variables with the same names override these.
# -----------------------------------------------------------------------------

# INPUT_SUFFIX makes the mapper reusable for the taxonomically narrowed read
# set without overwriting the existing whole-genus BAMs. For the new ORSC run:
#   BESTHIT_DIR=.../taxonomic_tiers INPUT_SUFFIX=.target_orsc.fastq.gz \
#   READSET_LABEL=target_orsc OUT_DIR=.../bam_irgsp_orsc ./map... submit all
BESTHIT_DIR="${BESTHIT_DIR:-/home/scratch/yinmt202607/gene/results/oryza_competitive_mapping/besthit}"
INPUT_SUFFIX="${INPUT_SUFFIX:-.besthit_oryza.fastq.gz}"
READSET_LABEL="${READSET_LABEL:-besthit_oryza}"
IRGSP_FA="${IRGSP_FA:-/home/scratch/yinmt202607/db/asian_rice_panel_index/irgsp.fa}"
OUT_DIR="${OUT_DIR:-/home/scratch/yinmt202607/gene/results/ecotype_pca/bam_irgsp}"
LOG_DIR="${OUT_DIR}/logs"
SUBMIT_DIR="${OUT_DIR}/submissions"

SLURM_PARTITION="${SLURM_PARTITION:-comp}"
SLURM_ACCOUNT="${SLURM_ACCOUNT:-}"

JOB_CPUS="${JOB_CPUS:-4}"
JOB_MEM_MB="${JOB_MEM_MB:-4000}"
JOB_TIME="${JOB_TIME:-00:30:00}"

SCRIPT_PATH="$(readlink -f "${BASH_SOURCE[0]}")"

usage() {
    cat <<'EOF'
Usage:
  map_besthit_to_irgsp.sh check
      Validate BESTHIT_DIR/IRGSP_FA paths, bwa+samtools on PATH (loading
      the samtools module if needed), and the SLURM partition/account.
      Submits nothing.

  map_besthit_to_irgsp.sh smoke [SAMPLE]
      Submit one job for one sample (first discovered if omitted).

  map_besthit_to_irgsp.sh submit SAMPLE [SAMPLE ...]
      Submit one full sbatch job per sample. Skips a sample if
      OUT_DIR/<sample>.finished already exists.

  map_besthit_to_irgsp.sh submit all
      Same, but auto-discovers every sample with a
      <sample><INPUT_SUFFIX> under BESTHIT_DIR.

  map_besthit_to_irgsp.sh local SAMPLE [SAMPLE ...] | local all
      Same worker as submit, but sequential in the foreground -- no
      sbatch, no queue wait. Aborts on the first failing sample.

  map_besthit_to_irgsp.sh run SAMPLE OUTDIR
      Internal: the actual per-sample worker (sbatch job body). Not
      meant to be called directly by a user.
EOF
}

report_error_reset() { trap report_error ERR; }

validate_readset_config() {
    [[ -n "$INPUT_SUFFIX" ]] || { echo "ERROR: INPUT_SUFFIX must not be empty" >&2; return 1; }
    [[ "$READSET_LABEL" =~ ^[A-Za-z0-9_.-]+$ ]] || {
        echo "ERROR: READSET_LABEL contains unsafe filename characters: $READSET_LABEL" >&2
        return 1
    }
}

timestamp() {
    date --iso-8601=seconds 2>/dev/null || date '+%Y-%m-%dT%H:%M:%S%z'
}

load_tools() {
    if command -v module >/dev/null 2>&1; then
        module load samtools/ >/dev/null 2>&1 || module load samtools >/dev/null 2>&1 || true
    fi
    command -v samtools >/dev/null 2>&1 || {
        echo "ERROR: samtools not found on PATH, even after 'module load samtools'." >&2
        echo "  Check 'module avail samtools' for the exact module name on this cluster." >&2
        exit 1
    }
    command -v bwa >/dev/null 2>&1 || {
        echo "ERROR: bwa not found on PATH (expected already on PATH via the base" >&2
        echo "  conda env on this cluster -- confirm the right env/module is active)." >&2
        exit 1
    }
}

sbatch_common_args() {
    SBATCH_COMMON=(--partition="$SLURM_PARTITION")
    if [[ -n "$SLURM_ACCOUNT" ]]; then
        SBATCH_COMMON+=(--account="$SLURM_ACCOUNT")
    fi
}

fastq_path() {
    printf '%s/%s%s\n' "$BESTHIT_DIR" "$1" "$INPUT_SUFFIX"
}

discover_samples() {
    validate_readset_config
    [[ -d "$BESTHIT_DIR" ]] || {
        echo "ERROR: BESTHIT_DIR not found: $BESTHIT_DIR" >&2
        return 1
    }
    shopt -s nullglob
    local files=("$BESTHIT_DIR"/*"$INPUT_SUFFIX")
    shopt -u nullglob
    SAMPLES=()
    local f base
    for f in "${files[@]}"; do
        base="$(basename "$f")"
        SAMPLES+=("${base%${INPUT_SUFFIX}}")
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
        --job-name=mapirgsp.preflight --cpus-per-task="$JOB_CPUS" \
        --mem="${JOB_MEM_MB}M" --time="$JOB_TIME" --wrap=true >/dev/null
}

run_check() {
    echo "[check] BESTHIT_DIR=$BESTHIT_DIR"
    echo "[check] INPUT_SUFFIX=$INPUT_SUFFIX READSET_LABEL=$READSET_LABEL"
    echo "[check] IRGSP_FA=$IRGSP_FA"
    echo "[check] OUT_DIR=$OUT_DIR"
    echo "[check] partition=$SLURM_PARTITION account=${SLURM_ACCOUNT:-<default>}"

    validate_readset_config

    [[ -f "$IRGSP_FA" ]] || { echo "ERROR: IRGSP_FA missing: $IRGSP_FA" >&2; exit 1; }

    discover_samples
    [[ "${#SAMPLES[@]}" -gt 0 ]] || {
        echo "ERROR: no *${INPUT_SUFFIX} under $BESTHIT_DIR" >&2
        exit 1
    }
    echo "[check] samples with a besthit FASTQ: ${#SAMPLES[@]}"
    printf '  %s\n' "${SAMPLES[@]}"

    load_tools
    validate_slurm_request
    echo "[check] PASS"
}

run_worker() {
    # Internal: invoked as the sbatch job body. Args: sample outdir
    [[ "$#" -ge 2 ]] || { usage >&2; exit 2; }
    local sample="$1" outdir="$2"
    validate_readset_config
    local fq
    fq="$(fastq_path "$sample")"
    [[ -f "$fq" ]] || { echo "ERROR: besthit FASTQ missing for $sample: $fq" >&2; exit 1; }

    load_tools
    mkdir -p "$outdir"

    echo "[run] sample=$sample started=$(timestamp)"
    echo "[run] fastq=$fq"
    echo "[run] outdir=$outdir"

    echo "[map] $sample: bwa aln + samse + filter (-F 0x904) + sort"
    bwa aln -l 1024 -n 0.01 -o 2 -t "$JOB_CPUS" "$IRGSP_FA" "$fq" > "$outdir/${sample}.sai"
    bwa samse "$IRGSP_FA" "$outdir/${sample}.sai" "$fq" \
        | samtools view -@ "$JOB_CPUS" -bh -F 0x904 - \
        | samtools sort -@ "$JOB_CPUS" -o "$outdir/${sample}.sorted.bam" -
    rm -f "$outdir/${sample}.sai"

    echo "[dedup] $sample: collate | fixmate -m | sort | markdup (flag only, not -r)"
    samtools collate -@ "$JOB_CPUS" -O "$outdir/${sample}.sorted.bam" \
        | samtools fixmate -@ "$JOB_CPUS" -m - - \
        | samtools sort -@ "$JOB_CPUS" -o "$outdir/${sample}.fixmate_sorted.bam" -
    local final_bam="$outdir/${sample}.${READSET_LABEL}.irgsp.bam"
    samtools markdup -@ "$JOB_CPUS" "$outdir/${sample}.fixmate_sorted.bam" "$final_bam"
    samtools index "$final_bam"
    rm -f "$outdir/${sample}.sorted.bam" "$outdir/${sample}.fixmate_sorted.bam"

    local n_mapped n_dup n_q30 n_q20
    n_mapped=$(samtools view -@ "$JOB_CPUS" -c -F 4 "$final_bam")
    n_dup=$(samtools view -@ "$JOB_CPUS" -c -f 1024 "$final_bam")
    n_q30=$(samtools view -@ "$JOB_CPUS" -c -F 1028 -q 30 "$final_bam")
    n_q20=$(samtools view -@ "$JOB_CPUS" -c -F 1028 -q 20 "$final_bam")
    echo "[done] $sample: mapped=$n_mapped duplicates_flagged=$n_dup mapq>=30_nondup=$n_q30 mapq>=20_nondup=$n_q20 (duplicates not removed -- pseudo_haploid_call.py filters them at pileup time)"

    touch "$outdir/${sample}.finished"
    echo "[run] sample=$sample finished=$(timestamp)"
}

run_smoke() {
    load_tools
    discover_samples
    local sample="${1:-}"
    if [[ -z "$sample" ]]; then
        [[ "${#SAMPLES[@]}" -gt 0 ]] || {
            echo "ERROR: no samples to pick a default from" >&2
            exit 1
        }
        sample="${SAMPLES[0]}"
    fi
    [[ -f "$(fastq_path "$sample")" ]] || {
        echo "ERROR: besthit FASTQ missing for $sample: $(fastq_path "$sample")" >&2
        exit 1
    }

    mkdir -p "$OUT_DIR" "$LOG_DIR" "$SUBMIT_DIR"
    validate_slurm_request
    sbatch_common_args

    local log="${LOG_DIR}/${sample}.smoke.%j.log"
    local job_id
    job_id="$(sbatch \
        --parsable \
        "${SBATCH_COMMON[@]}" \
        --job-name="mapirgsp.smoke.${sample}" \
        --cpus-per-task="$JOB_CPUS" \
        --mem="${JOB_MEM_MB}M" \
        --time="$JOB_TIME" \
        --output="$log" \
        --error="$log" \
        "$SCRIPT_PATH" run "$sample" "$OUT_DIR")"

    echo "[smoke] submitted job_id=$job_id sample=$sample"
    echo "[smoke] watch: squeue -j $job_id"
    echo "[smoke] log:   tail -f ${log//%j/$job_id}"
}

run_submit() {
    [[ "$#" -ge 1 ]] || { usage >&2; exit 2; }
    load_tools
    mkdir -p "$OUT_DIR" "$LOG_DIR" "$SUBMIT_DIR"
    validate_slurm_request
    sbatch_common_args

    local manifest="${SUBMIT_DIR}/submit.$(date +%Y%m%dT%H%M%S).tsv"
    printf 'job_id\tsample\n' > "$manifest"

    local sample fq log job_id
    for sample in "$@"; do
        fq="$(fastq_path "$sample")"
        [[ -f "$fq" ]] || { echo "ERROR: besthit FASTQ missing for $sample: $fq" >&2; exit 1; }

        if [[ -f "${OUT_DIR}/${sample}.finished" ]]; then
            echo "[submit] $sample already finished, skipping"
            continue
        fi

        log="${LOG_DIR}/${sample}.mapirgsp.%j.log"
        job_id="$(sbatch \
            --parsable \
            "${SBATCH_COMMON[@]}" \
            --job-name="mapirgsp.${sample}" \
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
    [[ "$#" -ge 1 ]] || { usage >&2; exit 2; }
    load_tools
    mkdir -p "$OUT_DIR" "$LOG_DIR"

    local sample fq log
    for sample in "$@"; do
        fq="$(fastq_path "$sample")"
        [[ -f "$fq" ]] || { echo "ERROR: besthit FASTQ missing for $sample: $fq" >&2; exit 1; }

        if [[ -f "${OUT_DIR}/${sample}.finished" ]]; then
            echo "[local] $sample already finished, skipping"
            continue
        fi

        log="${LOG_DIR}/${sample}.mapirgsp.local.log"
        echo "[local] running sample=$sample (foreground, no sbatch) -> $log"
        run_worker "$sample" "$OUT_DIR" 2>&1 | tee "$log"
    done
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
                echo "ERROR: no *${INPUT_SUFFIX} under $BESTHIT_DIR" >&2
                exit 1
            }
            echo "[submit] all: ${#SAMPLES[@]} samples with a besthit FASTQ"
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
                echo "ERROR: no *${INPUT_SUFFIX} under $BESTHIT_DIR" >&2
                exit 1
            }
            echo "[local] all: ${#SAMPLES[@]} samples with a besthit FASTQ"
            run_local "${SAMPLES[@]}"
        else
            run_local "$@"
        fi
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
