#!/usr/bin/env bash
# Batch classification-power validation for a given ancient sample's SNP
# mask, reusing the EXACT same chain already validated in
# docs/ECOTYPE_PCA_PHASE1_COMMANDS.md section 6 (single indica held-out
# LOO test) -- this script just loops that chain over many held-out
# individuals across several known population labels, instead of one.
#
# For each (label, held-out individual) pair: mask that individual to
# MASK_FROM's covered-site pattern (simulate_leaveoneout_projection.py),
# subset the panel to those same rows (build_sample_panel_subset.py
# --mask-from), relabel the held-out row so it can't help build its own
# axis, merge the simulated individual back in as a projected-only
# "Ancient"-labeled column (merge_ancient_into_panel.py), run smartpca
# with REFERENCE_LABELS_FILE restricting axis-building to the given
# domesticated labels, then rank its nearest population
# (summarize_projection_distances.py). No new classification logic is
# introduced anywhere in this chain -- it is the same nearest-centroid
# ranking already used for real ancient samples, just run on individuals
# whose true label is already known so the ranking's accuracy can be
# checked.
#
# Usage:
#   run_masked_loo_validation.sh MASK_FROM PANEL_SNP PANEL_GENO PANEL_IND \
#     OUT_DIR REFERENCE_LABELS_FILE N_PER_LABEL SEED LABEL [LABEL...]
#
# MASK_FROM: a pseudo_haploid_call.py --out file (e.g. the real ancient
#   sample's own pseudohap.txt) whose covered/missing SNP pattern every
#   held-out test individual is masked to.
# LABEL...: population labels (as they appear in PANEL_IND column 3) to
#   draw held-out test individuals from, e.g. aromatic japonica indica aus.
#
# Skips (does not re-run) any (label, held-out id) pair whose
# <label>.<id>.nearest.tsv already exists, so this can be re-submitted
# safely after a partial failure. Writes OUT_DIR/manifest.tsv, one row
# per (true_label, held_out_id, output_prefix) -- the input to
# build_confusion_matrix.py.

set -euo pipefail

if [[ $# -lt 9 ]]; then
    echo "usage: $0 MASK_FROM PANEL_SNP PANEL_GENO PANEL_IND OUT_DIR REFERENCE_LABELS_FILE N_PER_LABEL SEED LABEL [LABEL...]" >&2
    exit 1
fi

MASK_FROM="$1"
PANEL_SNP="$2"
PANEL_GENO="$3"
PANEL_IND="$4"
OUT_DIR="$5"
REFERENCE_LABELS_FILE="$6"
N_PER_LABEL="$7"
SEED="$8"
shift 8
LABELS=("$@")

for f in "$MASK_FROM" "$PANEL_SNP" "$PANEL_GENO" "$PANEL_IND" "$REFERENCE_LABELS_FILE"; do
    [[ -s "$f" ]] || { echo "ERROR: missing or empty input file: $f" >&2; exit 1; }
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "$OUT_DIR"

MANIFEST="${OUT_DIR}/manifest.tsv"
if [[ ! -s "$MANIFEST" ]]; then
    printf 'true_label\theld_out_id\tprefix\n' > "$MANIFEST"
fi

PROCESSED_PREFIXES=()

for LABEL in "${LABELS[@]}"; do
    N_AVAILABLE=$(awk -v lab="$LABEL" '$3==lab' "$PANEL_IND" | wc -l)
    if [[ "$N_AVAILABLE" -eq 0 ]]; then
        echo "ERROR: no individuals with label '${LABEL}' found in ${PANEL_IND}" >&2
        exit 1
    fi

    SELECTED_IDS=$(awk -v lab="$LABEL" '$3==lab{print $1}' "$PANEL_IND" \
        | python3 -c "
import random, sys
ids = sys.stdin.read().split()
random.Random(${SEED}).shuffle(ids)
print('\n'.join(ids[:${N_PER_LABEL}]))
")

    for HELD_OUT in $SELECTED_IDS; do
        PREFIX="${OUT_DIR}/${LABEL}.${HELD_OUT}"

        if [[ -s "${PREFIX}.nearest.tsv" ]]; then
            echo "[skip] ${PREFIX}.nearest.tsv already exists"
            continue
        fi

        echo "[run] true_label=${LABEL} held_out=${HELD_OUT}"

        python3 "${SCRIPT_DIR}/simulate_leaveoneout_projection.py" \
            --panel-geno "$PANEL_GENO" --panel-ind "$PANEL_IND" \
            --held-out-sample "$HELD_OUT" --mask-from "$MASK_FROM" --seed "$SEED" \
            --out "${PREFIX}.simulated.calls.txt" \
            --report "${PREFIX}.simulated.report.tsv"

        python3 "${SCRIPT_DIR}/build_sample_panel_subset.py" \
            --panel-snp "$PANEL_SNP" --panel-geno "$PANEL_GENO" \
            --ancient-calls "${PREFIX}.simulated.calls.txt" --mask-from "$MASK_FROM" \
            --out-panel-snp "${PREFIX}.subset.snp" \
            --out-panel-geno "${PREFIX}.subset.eigenstratgeno" \
            --out-ancient-calls "${PREFIX}.subset.calls.txt" \
            --report "${PREFIX}.subset.report.tsv"

        awk -v s="$HELD_OUT" '{if($1==s){$3="LOO_HELDOUT_EXCLUDED"}; print}' "$PANEL_IND" \
            > "${PREFIX}.loo_test.ind"

        python3 "${SCRIPT_DIR}/merge_ancient_into_panel.py" \
            --panel-geno "${PREFIX}.subset.eigenstratgeno" \
            --panel-ind "${PREFIX}.loo_test.ind" \
            --calls "${HELD_OUT}_LOOSIM=${PREFIX}.subset.calls.txt" \
            --ancient-poplabel Ancient \
            --out-geno "${PREFIX}.merged.eigenstratgeno" \
            --out-ind "${PREFIX}.merged.ind"

        sort -u "$REFERENCE_LABELS_FILE" > "${PREFIX}.poplistname.txt"

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

        python3 "${SCRIPT_DIR}/summarize_projection_distances.py" \
            --evec "${PREFIX}.evec" --sample "${HELD_OUT}_LOOSIM" \
            --out "${PREFIX}.nearest.tsv"

        printf '%s\t%s\t%s\n' "$LABEL" "$HELD_OUT" "$PREFIX" >> "$MANIFEST"
        PROCESSED_PREFIXES+=("$PREFIX")
    done
done

# Sanity check: row selection is driven only by MASK_FROM's covered
# positions, not by which individual is held out -- every *.subset.snp
# produced in this run should therefore be byte-identical. A mismatch
# means something is wrong with the row-alignment assumption itself.
if [[ "${#PROCESSED_PREFIXES[@]}" -gt 1 ]]; then
    DISTINCT_SNP_MD5=$(md5sum "${PROCESSED_PREFIXES[@]/%/.subset.snp}" | awk '{print $1}' | sort -u | wc -l)
    if [[ "$DISTINCT_SNP_MD5" -gt 1 ]]; then
        echo "WARNING: this run's *.subset.snp files are not all identical -- row alignment assumption may be violated, check manually" >&2
    fi
fi

echo "[done] manifest: $MANIFEST"
