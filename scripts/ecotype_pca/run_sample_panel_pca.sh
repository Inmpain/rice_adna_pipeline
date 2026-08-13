#!/usr/bin/env bash
# Orchestrates ONE sample x ONE panel through the full sample-specific
# subset PCA chain (docs/ECOTYPE_PCA_EXECUTION_PLAN.md section 5):
#   pseudo_haploid_call.py -> build_sample_panel_subset.py ->
#   merge_ancient_into_panel.py -> smartpca -lsqproject
#
# Every ancient sample gets its OWN subset panel (shrunk to just the
# SNPs this one sample's reads cover), so this is run once per
# sample-panel pair, not once for a batch of samples merged together --
# section 5's PCA coordinates are explicitly NOT comparable across
# samples run this way.
#
# Usage:
#   run_sample_panel_pca.sh SAMPLE PANEL_LABEL BAM PANEL_SNP PANEL_GENO PANEL_IND OUT_DIR [TV|ALL] [REFERENCE_LABELS_FILE]
#
# PANEL_GENO/PANEL_IND must already be the real-population-labeled,
# UNK-dropped versions (build_*_population_labels.py + filter_panel_by_label.py
# output), not the raw panel files -- this script does not do that step.
#
# REFERENCE_LABELS_FILE (optional, 9th arg): a file listing (one per
# line) exactly which population labels are allowed to build the PCA
# axes -- passed straight through as smartpca's poplistname. Any label
# in PANEL_IND that is NOT in this file still gets merged in and
# -lsqproject'd onto the resulting axes, just like the ancient sample,
# it just doesn't influence where those axes are.
#
# WHY THIS EXISTS (2026-08-13, found via GPT review): without it, this
# script's old behavior put EVERY non-Ancient label into poplistname --
# on the Civáň panel that means the 456 O._rufipogon wild samples (plus
# 4 singleton wild outgroups) helped define the axes the ancient sample
# gets projected onto. Civáň et al. 2019's own method is the opposite:
# build axes from the 595 domesticated accessions only, project wild
# rice onto them (justified in the paper by wild rice's much higher
# missingness -- letting it into axis-building would let coverage
# artifacts, not just ancestry, shape the PCs). Omit this argument to
# keep the old include-everything behavior (fine for 29M_3k, which is
# 100% cultivated already; 6.7M_720's correct reference set is a
# separate, not-yet-decided design question -- do not assume
# domesticated-only is the right call there too).
#
# Skips (does not re-run) any sample-panel-track combination whose
# .evec already exists and is non-empty, so a batch loop can be
# re-submitted safely after a partial failure.

set -euo pipefail

if [[ $# -lt 7 ]]; then
    echo "usage: $0 SAMPLE PANEL_LABEL BAM PANEL_SNP PANEL_GENO PANEL_IND OUT_DIR [TV|ALL] [REFERENCE_LABELS_FILE]" >&2
    exit 1
fi

SAMPLE="$1"
PANEL_LABEL="$2"
BAM="$3"
PANEL_SNP="$4"
PANEL_GENO="$5"
PANEL_IND="$6"
OUT_DIR="$7"
TRACK="${8:-TV}"
REFERENCE_LABELS_FILE="${9:-}"

if [[ "$TRACK" != "TV" && "$TRACK" != "ALL" ]]; then
    echo "ERROR: track must be TV or ALL, got: $TRACK" >&2
    exit 1
fi

if [[ -n "$REFERENCE_LABELS_FILE" ]]; then
    [[ -s "$REFERENCE_LABELS_FILE" ]] || { echo "ERROR: --reference-labels-file given but missing or empty: $REFERENCE_LABELS_FILE" >&2; exit 1; }
fi

for f in "$BAM" "$PANEL_SNP" "$PANEL_GENO" "$PANEL_IND"; do
    [[ -s "$f" ]] || { echo "ERROR: missing or empty input file: $f" >&2; exit 1; }
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "$OUT_DIR"

PREFIX="${OUT_DIR}/${SAMPLE}.${PANEL_LABEL}.${TRACK}"

if [[ -s "${PREFIX}.evec" ]]; then
    echo "[skip] ${PREFIX}.evec already exists"
    exit 0
fi

CALL_ARGS=()
if [[ "$TRACK" == "ALL" ]]; then
    CALL_ARGS+=(--no-transversions-only)
fi

echo "[1/4] pseudo_haploid_call.py (${TRACK})"
python3 "${SCRIPT_DIR}/pseudo_haploid_call.py" \
    --bam "$BAM" \
    --panel-snp "$PANEL_SNP" \
    --out "${PREFIX}.pseudohap.txt" \
    --report "${PREFIX}.pseudohap.report.tsv" \
    "${CALL_ARGS[@]}"

echo "[2/4] build_sample_panel_subset.py"
python3 "${SCRIPT_DIR}/build_sample_panel_subset.py" \
    --panel-snp "$PANEL_SNP" \
    --panel-geno "$PANEL_GENO" \
    --ancient-calls "${PREFIX}.pseudohap.txt" \
    --out-panel-snp "${PREFIX}.subset.snp" \
    --out-panel-geno "${PREFIX}.subset.eigenstratgeno" \
    --out-ancient-calls "${PREFIX}.subset.calls.txt" \
    --report "${PREFIX}.subset.report.tsv"

echo "[3/4] merge_ancient_into_panel.py"
python3 "${SCRIPT_DIR}/merge_ancient_into_panel.py" \
    --panel-geno "${PREFIX}.subset.eigenstratgeno" \
    --panel-ind "$PANEL_IND" \
    --calls "${SAMPLE}=${PREFIX}.subset.calls.txt" \
    --ancient-poplabel Ancient \
    --out-geno "${PREFIX}.merged.eigenstratgeno" \
    --out-ind "${PREFIX}.merged.ind"

echo "[4/4] smartpca"
if [[ -n "$REFERENCE_LABELS_FILE" ]]; then
    sort -u "$REFERENCE_LABELS_FILE" > "${PREFIX}.poplistname.txt"
    PRESENT_LABELS="$(awk '{print $3}' "${PREFIX}.merged.ind" | sort -u)"
    MISSING_LABELS="$(comm -23 "${PREFIX}.poplistname.txt" <(echo "$PRESENT_LABELS"))"
    if [[ -n "$MISSING_LABELS" ]]; then
        echo "ERROR: --reference-labels-file lists label(s) with zero individuals in ${PREFIX}.merged.ind: $(echo "$MISSING_LABELS" | tr '\n' ' ')" >&2
        exit 1
    fi
else
    awk '{print $3}' "${PREFIX}.merged.ind" | sort -u | grep -v "^Ancient$" > "${PREFIX}.poplistname.txt"
fi

cat > "${PREFIX}.par" << PAREOF
genotypename:    ${PREFIX}.merged.eigenstratgeno
snpname:         ${PREFIX}.subset.snp
indivname:       ${PREFIX}.merged.ind
evecoutname:     ${PREFIX}.evec
evaloutname:     ${PREFIX}.eval
poplistname:     ${PREFIX}.poplistname.txt
lsqproject:      YES
numoutevec:      10
numoutlieriter:  0
numchrom:        12
numthreads:      2
PAREOF

smartpca -p "${PREFIX}.par" > "${PREFIX}.smartpca.log" 2>&1

if [[ ! -s "${PREFIX}.evec" ]]; then
    echo "ERROR: smartpca did not produce ${PREFIX}.evec -- see ${PREFIX}.smartpca.log" >&2
    exit 1
fi

echo "[done] ${PREFIX}.evec"
