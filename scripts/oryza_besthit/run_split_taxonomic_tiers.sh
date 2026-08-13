#!/usr/bin/env bash
# Split existing genus-wide best-hit outputs without rerunning competitive mapping.
set -Eeuo pipefail

BESTHIT_DIR="${BESTHIT_DIR:-/home/scratch/yinmt202607/gene/results/oryza_competitive_mapping/besthit}"
OUT_DIR="${OUT_DIR:-/home/scratch/yinmt202607/gene/results/oryza_competitive_mapping/taxonomic_tiers}"
TARGET_TAXIDS="${TARGET_TAXIDS:-4529,4530,4536}"
TARGET_LABEL="${TARGET_LABEL:-target_orsc}"

SCRIPT_PATH="$(readlink -f "${BASH_SOURCE[0]}")"
SCRIPT_DIR="$(dirname "$SCRIPT_PATH")"
PY_SCRIPT="${PY_SCRIPT:-${SCRIPT_DIR}/split_besthit_taxonomic_tiers.py}"

usage() {
    cat <<'EOF'
Usage:
  run_split_taxonomic_tiers.sh check
  run_split_taxonomic_tiers.sh run SAMPLE [SAMPLE ...]
  run_split_taxonomic_tiers.sh run all
  run_split_taxonomic_tiers.sh merge

This stage is small and I/O-bound, so it intentionally runs locally rather
than submitting 16 tiny SLURM jobs. Existing completed outputs are skipped.

Environment overrides:
  BESTHIT_DIR, OUT_DIR, TARGET_TAXIDS, TARGET_LABEL, PY_SCRIPT
EOF
}

discover_samples() {
    [[ -d "$BESTHIT_DIR" ]] || {
        echo "ERROR: BESTHIT_DIR not found: $BESTHIT_DIR" >&2
        return 1
    }
    shopt -s nullglob
    local files=("$BESTHIT_DIR"/*.besthit_oryza.fastq.gz)
    shopt -u nullglob
    SAMPLES=()
    local file base
    for file in "${files[@]}"; do
        base="$(basename "$file")"
        SAMPLES+=("${base%.besthit_oryza.fastq.gz}")
    done
}

validate_sample() {
    local sample="$1"
    [[ -f "$BESTHIT_DIR/${sample}.besthit_oryza.fastq.gz" ]] || {
        echo "ERROR: missing best-hit FASTQ for $sample" >&2
        return 1
    }
    [[ -f "$BESTHIT_DIR/${sample}.oryza_filter.decisions.tsv.gz" ]] || {
        echo "ERROR: missing decisions table for $sample" >&2
        return 1
    }
}

run_check() {
    command -v python3 >/dev/null || { echo "ERROR: python3 not found" >&2; exit 1; }
    [[ -f "$PY_SCRIPT" ]] || { echo "ERROR: PY_SCRIPT not found: $PY_SCRIPT" >&2; exit 1; }
    discover_samples
    [[ "${#SAMPLES[@]}" -gt 0 ]] || {
        echo "ERROR: no *.besthit_oryza.fastq.gz under $BESTHIT_DIR" >&2
        exit 1
    }
    local sample
    for sample in "${SAMPLES[@]}"; do
        validate_sample "$sample"
    done
    echo "[check] PASS samples=${#SAMPLES[@]} target_taxids=$TARGET_TAXIDS"
    echo "[check] output=$OUT_DIR"
}

run_samples() {
    mkdir -p "$OUT_DIR"
    local sample target_out
    for sample in "$@"; do
        validate_sample "$sample"
        target_out="$OUT_DIR/${sample}.${TARGET_LABEL}.fastq.gz"
        if [[ -f "$target_out" && -f "$OUT_DIR/${sample}.taxonomic_tiers.summary.tsv" ]]; then
            echo "[run] $sample already complete, skipping"
            continue
        fi
        echo "[run] splitting $sample"
        python3 "$PY_SCRIPT" \
            --sample "$sample" \
            --besthit-fastq "$BESTHIT_DIR/${sample}.besthit_oryza.fastq.gz" \
            --decisions "$BESTHIT_DIR/${sample}.oryza_filter.decisions.tsv.gz" \
            --outdir "$OUT_DIR" \
            --target-taxids "$TARGET_TAXIDS" \
            --target-label "$TARGET_LABEL"
    done
}

run_merge() {
    shopt -s nullglob
    local summaries=("$OUT_DIR"/*.taxonomic_tiers.summary.tsv)
    local species=("$OUT_DIR"/*.taxonomic_tiers.by_species.tsv)
    shopt -u nullglob
    [[ "${#summaries[@]}" -gt 0 ]] || { echo "ERROR: no summaries under $OUT_DIR" >&2; exit 1; }
    [[ "${#species[@]}" -eq "${#summaries[@]}" ]] || {
        echo "ERROR: summary/by-species file counts differ" >&2
        exit 1
    }

    local summary_out="$OUT_DIR/taxonomic_tiers_summary.tsv"
    local species_out="$OUT_DIR/taxonomic_tiers_by_species.tsv"
    local summary_tmp="${summary_out}.tmp.$$"
    local species_tmp="${species_out}.tmp.$$"
    trap 'rm -f "$summary_tmp" "$species_tmp"' EXIT

    head -n 1 "${summaries[0]}" > "$summary_tmp"
    head -n 1 "${species[0]}" > "$species_tmp"
    local file
    for file in "${summaries[@]}"; do tail -n +2 "$file" >> "$summary_tmp"; done
    for file in "${species[@]}"; do tail -n +2 "$file" >> "$species_tmp"; done
    mv "$summary_tmp" "$summary_out"
    mv "$species_tmp" "$species_out"
    trap - EXIT
    echo "[merge] ${#summaries[@]} samples -> $summary_out"
    echo "[merge] species detail -> $species_out"
}

case "${1:-}" in
    check)
        run_check
        ;;
    run)
        shift
        [[ "$#" -gt 0 ]] || { usage >&2; exit 2; }
        if [[ "${1:-}" == "all" ]]; then
            discover_samples
            run_samples "${SAMPLES[@]}"
        else
            run_samples "$@"
        fi
        ;;
    merge)
        run_merge
        ;;
    *)
        usage
        exit 2
        ;;
esac
