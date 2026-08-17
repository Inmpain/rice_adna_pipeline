#!/usr/bin/env bash
# Auxiliary/cross-check shared-matrix projection for Panel A (3K) and Panel B
# (720) -- per docs/ECOTYPE_PCA_EXECUTION_PLAN.md section 5's "辅助分析"
# (auxiliary, NOT the primary analysis): a single MAF-primary, lightly-
# LD-pruned modern reference marker set, all 16 ancient samples lsqproject'd
# onto the SAME axes so they can be plotted together on one PC1-10 grid.
# The primary, per-sample analysis is v1's scripts/ecotype_pca/
# run_sample_panel_pca.sh (each sample against its own real full-panel
# coverage, coordinates not comparable across samples) -- this script does
# not replace that, it is the secondary cross-check.
#
# Deliberately simple, unlike 51_civan_maf_ld_and_private_axis.sh: no
# ancient-coverage-first restriction, no TV track, no private axis here.
# Standard 07 "fixed"/primary route (geno -> MAF -> LD, values from config)
# is used as-is; whatever marker count survives is what gets used.
set -euo pipefail
: "${SLURM_JOB_ID:?stage 61 must run in SLURM}"
: "${RICE_PCA_REPO_ROOT:?}"; : "${RICE_PCA_CONFIG:?}"; : "${RICE_PCA_ATTEMPT_DIR:?}"
: "${PANEL_LETTER:?A or B}"
: "${PANEL_ANCIENT_SAMPLES:?space-separated list of ancient sample IDs}"
cd "$RICE_PCA_REPO_ROOT"

case "$PANEL_LETTER" in
  A) PANEL_KEY=panel_A_3k; PANEL_LABEL=3k ;;
  B) PANEL_KEY=panel_B_720; PANEL_LABEL=720 ;;
  *) echo "FATAL: PANEL_LETTER must be A or B, got '$PANEL_LETTER'" >&2; exit 1 ;;
esac

eval "$(python3 - "$RICE_PCA_CONFIG" "$PANEL_KEY" <<'PY'
import shlex,sys,yaml
c=yaml.safe_load(open(sys.argv[1])); key=sys.argv[2]; p=c['inputs'][key]
geno_ext = p.get('geno_ext', '.eigenstratgeno')
for k,v in {
  'PANEL_DIR': p['dir'],
  'PANEL_PREFIX': p['prefix'],
  'PANEL_SNP': f"{p['dir']}/{p['prefix']}.snp",
  'PANEL_GENO': f"{p['dir']}/{p['prefix']}{p['filtered_suffix']}{geno_ext}",
  'PANEL_IND': f"{p['dir']}/{p['prefix']}{p['filtered_suffix']}.ind",
  'BAMDIR': c['inputs']['ancient_bam_dir'],
}.items(): print(f'{k}={shlex.quote(v)}')
PY
)"

OUT="$RICE_PCA_ATTEMPT_DIR"
mkdir -p "$OUT/plink_input_staging" "$OUT/plink" "$OUT/reference_sets" "$OUT/maf_ld" "$OUT/SHARED/calls"

echo "=== Step 1: stage $PANEL_LABEL raw .snp + filtered .ind/geno under one shared prefix ==="
# Same split-naming workaround as 51_civan_maf_ld_and_private_axis.sh: 02
# requires PREFIX.ind/.snp/.eigenstratgeno in one dir sharing one prefix;
# symlink rather than touch the real files (Panel B's real geno file is
# named .geno, not .eigenstratgeno -- handled by geno_ext above, the
# symlink target name is always .eigenstratgeno regardless).
STAGE_DIR="$OUT/plink_input_staging"
ln -sf "$PANEL_SNP" "$STAGE_DIR/${PANEL_PREFIX}.snp"
ln -sf "$PANEL_IND" "$STAGE_DIR/${PANEL_PREFIX}.ind"
ln -sf "$PANEL_GENO" "$STAGE_DIR/${PANEL_PREFIX}.eigenstratgeno"

bash scripts/ecotype_pca_v2/02_convert_eigenstrat_for_plink.sh \
  --dir "$STAGE_DIR" --prefix "$PANEL_PREFIX" --out-dir "$OUT/plink"
PANEL_BFILE="$OUT/plink/${PANEL_PREFIX}.plink"

echo "=== Step 1b: build panel $PANEL_LETTER reference/axis-builder keep-list (real FID from the .fam just produced) ==="
python3 scripts/ecotype_pca_v2/06_build_reference_sample_set.py \
  --config "$RICE_PCA_CONFIG" --panel "$PANEL_LETTER" --label "$PANEL_LABEL" \
  --ind-file "$PANEL_IND" --fam-file "${PANEL_BFILE}.fam" \
  --out-dir "$OUT/reference_sets"
REFERENCE_KEEP="$OUT/reference_sets/${PANEL_LABEL}.reference_samples.keep"

echo "=== MAF/LD marker selection (07, panel=$PANEL_LETTER, sensitivity=primary, ALL track) ==="
bash scripts/ecotype_pca_v2/07_make_fixed_markers.sh \
  --config "$RICE_PCA_CONFIG" --panel "$PANEL_LETTER" --sensitivity primary \
  --library-type pooled_mixed --track ALL \
  --bfile "$PANEL_BFILE" --keep "$REFERENCE_KEEP" \
  --label "$PANEL_LABEL" --out-dir "$OUT/maf_ld"
FIXED_SNPLIST="$OUT/maf_ld/${PANEL_LABEL}.pooled_mixed.ALL.primary.fixed.snplist"
echo "marker count: $(wc -l < "$FIXED_SNPLIST")"

echo "=== shared-matrix fixed reference on the MAF/LD-cleaned marker set (09) ==="
SHARED_DIR="$OUT/SHARED"
python3 scripts/ecotype_pca_v2/09_export_fixed_reference_eigenstrat.py \
  --panel "$PANEL_LETTER" --library-type pooled_mixed --track ALL \
  --panel-snp "$PANEL_SNP" --panel-geno "$PANEL_GENO" --panel-ind "$PANEL_IND" \
  --fixed-snplist "$FIXED_SNPLIST" --reference-keep "$REFERENCE_KEEP" \
  --label "$PANEL_LABEL" --out-dir "$SHARED_DIR"
REFERENCE_PREFIX="$SHARED_DIR/${PANEL_LABEL}.pooled_mixed.ALL.fixed_reference"

SUMMARY="$OUT/stage61_summary.tsv"
printf 'sample\ttechnical_execution\ttechnical_note\n' > "$SUMMARY"

CALLS_ARGS=()
for SAMPLE in $PANEL_ANCIENT_SAMPLES; do
  BAM="$BAMDIR/$SAMPLE.besthit_oryza.irgsp.bam"
  if [[ ! -s "$BAM" ]]; then
    echo "WARNING: $SAMPLE: BAM missing or empty, technical_execution=FAIL" >&2
    printf '%s\tFAIL\tBAM missing or empty\n' "$SAMPLE" >> "$SUMMARY"
    continue
  fi
  set +e
  python3 scripts/ecotype_pca_v2/10_call_ancient_fixed_markers.py \
    --config "$RICE_PCA_CONFIG" --bam "$BAM" --fixed-snp "$REFERENCE_PREFIX.snp" \
    --sample "$SAMPLE" --panel "$PANEL_LETTER" --library-type pooled_mixed --track ALL \
    --out-dir "$SHARED_DIR/calls" 2> "$SHARED_DIR/calls/$SAMPLE.stderr.log"
  RC=$?
  set -e
  if [[ $RC -eq 0 ]]; then
    printf '%s\tPASS\t\n' "$SAMPLE" >> "$SUMMARY"
    CALLS_ARGS+=("$SAMPLE=$SHARED_DIR/calls/$SAMPLE.$PANEL_LETTER.pooled_mixed.ALL.calls.txt")
  else
    NOTE=$(tail -1 "$SHARED_DIR/calls/$SAMPLE.stderr.log" 2>/dev/null || echo "unknown error")
    echo "WARNING: $SAMPLE: calling failed (exit $RC): $NOTE" >&2
    printf '%s\tFAIL\t%s\n' "$SAMPLE" "$NOTE" >> "$SUMMARY"
  fi
done
if [[ ${#CALLS_ARGS[@]} -eq 0 ]]; then
  echo "FATAL: no sample completed calling against the MAF/LD-cleaned marker set" >&2
  exit 3
fi

python3 scripts/ecotype_pca_v2/11_build_ancient_callability.py \
  --config "$RICE_PCA_CONFIG" --fixed-snp "$REFERENCE_PREFIX.snp" \
  $(for a in "${CALLS_ARGS[@]}"; do printf -- '--calls %s ' "$a"; done) \
  --panel "$PANEL_LETTER" --library-type pooled_mixed --track ALL \
  --out "$SHARED_DIR/${PANEL_LABEL}.callability.tsv"

python3 scripts/ecotype_pca_v2/13_merge_ancients_fixed_panel.py \
  --reference-geno "$REFERENCE_PREFIX.eigenstratgeno" --reference-ind "$REFERENCE_PREFIX.ind" \
  --fixed-snp "$REFERENCE_PREFIX.snp" \
  $(for a in "${CALLS_ARGS[@]}"; do printf -- '--calls %s ' "$a"; done) \
  --ancient-poplabel Ancient --label "$PANEL_LABEL" --out-dir "$SHARED_DIR"

bash scripts/ecotype_pca_v2/14_run_fixed_smartpca.sh \
  --config "$RICE_PCA_CONFIG" \
  --geno "$SHARED_DIR/${PANEL_LABEL}.merged.eigenstratgeno" --snp "$REFERENCE_PREFIX.snp" \
  --ind "$SHARED_DIR/${PANEL_LABEL}.merged.ind" --poplist "$REFERENCE_PREFIX.poplistname" \
  --label "${PANEL_LABEL}.pca" --out-dir "$SHARED_DIR"

EXPECTED_N=$(wc -l < "$SHARED_DIR/${PANEL_LABEL}.merged.ind")
python3 scripts/ecotype_pca_v2/15_pca_qc.py \
  --evec "$SHARED_DIR/${PANEL_LABEL}.pca.evec" --ind "$SHARED_DIR/${PANEL_LABEL}.merged.ind" \
  --expected-n "$EXPECTED_N" --out "$SHARED_DIR/${PANEL_LABEL}.pca_qc.tsv"

python3 scripts/ecotype_pca_v2/16_projection_summary.py \
  --panel "$PANEL_LETTER" --evec "$SHARED_DIR/${PANEL_LABEL}.pca.evec" \
  --out "$SHARED_DIR/${PANEL_LABEL}.projection_summary.tsv"

REPORT_ARGS=()
for a in "${CALLS_ARGS[@]}"; do
  SAMPLE="${a%%=*}"
  REPORT_ARGS+=("$SAMPLE=$SHARED_DIR/calls/$SAMPLE.$PANEL_LETTER.pooled_mixed.ALL.call_report.tsv")
done
python3 scripts/ecotype_pca_v2/22_classify_scientific_projection.py \
  $(for a in "${REPORT_ARGS[@]}"; do printf -- '--call-report %s ' "$a"; done) \
  --out "$SHARED_DIR/${PANEL_LABEL}.scientific_projection.tsv"

echo "=== plotting all PC1-PC10 pairs (26) ==="
bash scripts/ecotype_pca_v2/26_plot_pc_pairs.sh \
  --evec-label "${PANEL_LABEL}_shared" --evec "$SHARED_DIR/${PANEL_LABEL}.pca.evec" \
  --title-prefix "Panel $PANEL_LETTER ($PANEL_LABEL) MAF-primary shared matrix (ALL track)" \
  --out-prefix "$SHARED_DIR/${PANEL_LABEL}.shared"

echo "PASS: 61_panel_maf_shared_projection ($PANEL_LABEL) complete"
echo "See $SUMMARY, $SHARED_DIR/${PANEL_LABEL}.scientific_projection.tsv, and $SHARED_DIR/*.png"
