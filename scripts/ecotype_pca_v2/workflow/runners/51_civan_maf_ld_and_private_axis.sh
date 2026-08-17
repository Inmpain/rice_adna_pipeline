#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?stage 51 must run in SLURM}"
: "${RICE_PCA_REPO_ROOT:?}"; : "${RICE_PCA_CONFIG:?}"; : "${RICE_PCA_ATTEMPT_DIR:?}"
: "${CIVAN_UNION_SITES:?path to ancient_union_sites.tsv from 19_survey_ancient_coverage.py}"
: "${CIVAN_UNION_SITES_TV:?path to ancient_union_sites.TV.tsv from 20_filter_coverage_sites_to_transversions.py}"
: "${CIVAN_REFERENCE_FASTA:?path to the IRGSP reference FASTA the ancient BAMs were mapped against}"
: "${CIVAN_ANCIENT_SAMPLES:?space-separated list of ancient sample IDs}"
cd "$RICE_PCA_REPO_ROOT"

eval "$(python3 - "$RICE_PCA_CONFIG" <<'PY'
import shlex,sys,yaml
c=yaml.safe_load(open(sys.argv[1])); p=c['inputs']['panel_C_civan']
for k,v in {
  'PANEL_DIR': p['dir'],
  'PANEL_PREFIX': p['prefix'],
  'PANEL_SNP': f"{p['dir']}/{p['prefix']}.snp",
  'PANEL_GENO': f"{p['dir']}/{p['prefix']}{p['filtered_suffix']}.eigenstratgeno",
  'PANEL_IND': f"{p['dir']}/{p['prefix']}{p['filtered_suffix']}.ind",
  'BAMDIR': c['inputs']['ancient_bam_dir'],
}.items(): print(f'{k}={shlex.quote(v)}')
PY
)"

OUT="$RICE_PCA_ATTEMPT_DIR"
mkdir -p "$OUT/plink_input_staging" "$OUT/plink" "$OUT/reference_sets" "$OUT/maf_ld" "$OUT/SHARED/TV" "$OUT/SHARED/ALL" "$OUT/PRIVATE"

echo "=== Step 0: REF-vs-FASTA validation (once, against the raw Civan panel) ==="
python3 scripts/ecotype_pca_v2/23_validate_snp_ref_against_fasta.py \
  --snp "$PANEL_SNP" --fasta "$CIVAN_REFERENCE_FASTA" \
  --out "$OUT/civan.ref_vs_fasta.report.tsv"

echo "=== Step 1: stage Civan's raw .snp + filtered .ind/.eigenstratgeno under one shared prefix ==="
# 02_convert_eigenstrat_for_plink.sh requires PREFIX.ind / PREFIX.snp /
# PREFIX.eigenstratgeno to share one prefix in one directory. The real
# config splits these across an unfiltered .snp and a .filtered.ind /
# .filtered.eigenstratgeno -- stage symlinks rather than editing 02 (which
# is shared with panels A/B and has the identical split-naming issue there).
STAGE_DIR="$OUT/plink_input_staging"
ln -sf "$PANEL_SNP" "$STAGE_DIR/${PANEL_PREFIX}.snp"
ln -sf "$PANEL_IND" "$STAGE_DIR/${PANEL_PREFIX}.ind"
ln -sf "$PANEL_GENO" "$STAGE_DIR/${PANEL_PREFIX}.eigenstratgeno"

bash scripts/ecotype_pca_v2/02_convert_eigenstrat_for_plink.sh \
  --dir "$STAGE_DIR" --prefix "$PANEL_PREFIX" --out-dir "$OUT/plink"
CIVAN_BFILE="$OUT/plink/${PANEL_PREFIX}.plink"

echo "=== Step 1b: build panel C reference/axis-builder keep-list (real FID from the .fam just produced) ==="
# convertf's familynames: NO assigns sequential-index FIDs in the real
# .fam, not the sample ID string -- 06_build_reference_sample_set.py must
# be given --fam-file to write a keep-list plink2 --keep can actually
# match, or every row silently fails to match and plink2 keeps 0 samples.
python3 scripts/ecotype_pca_v2/06_build_reference_sample_set.py \
  --config "$RICE_PCA_CONFIG" --panel C --label civan \
  --ind-file "$PANEL_IND" --fam-file "${CIVAN_BFILE}.fam" \
  --out-dir "$OUT/reference_sets"
CIVAN_REFERENCE_KEEP="$OUT/reference_sets/civan.reference_samples.keep"

SUMMARY="$OUT/stage51_summary.tsv"
printf 'section\ttrack\tsample\ttechnical_execution\ttechnical_note\tcallable_n\n' > "$SUMMARY"

for TRACK in ALL TV; do
  echo "=== Track $TRACK: MAF/LD marker selection (07) ==="
  bash scripts/ecotype_pca_v2/07_make_fixed_markers.sh \
    --config "$RICE_PCA_CONFIG" --panel C --sensitivity primary \
    --library-type pooled_mixed --track "$TRACK" \
    --bfile "$CIVAN_BFILE" --keep "$CIVAN_REFERENCE_KEEP" \
    --label civan --out-dir "$OUT/maf_ld"
  MAF_LD_SNPLIST="$OUT/maf_ld/civan.pooled_mixed.${TRACK}.primary.fixed.snplist"

  echo "=== Track $TRACK: intersect MAF/LD-clean IDs with ancient-coverage IDs ==="
  case "$TRACK" in
    ALL) COVERAGE_SITES="$CIVAN_UNION_SITES" ;;
    TV)  COVERAGE_SITES="$CIVAN_UNION_SITES_TV" ;;
  esac
  TRACK_DIR="$OUT/SHARED/$TRACK"
  mkdir -p "$TRACK_DIR/calls"
  COVERAGE_SNPLIST="$TRACK_DIR/coverage.snplist.txt"
  python3 scripts/ecotype_pca_v2/21_extract_fixed_snplist.py \
    --sites-tsv "$COVERAGE_SITES" --out "$COVERAGE_SNPLIST"

  CLEAN_SNPLIST="$TRACK_DIR/civan.maf_ld_and_coverage.fixed.snplist"
  python3 scripts/ecotype_pca_v2/25_intersect_snplists.py \
    --snplist "$MAF_LD_SNPLIST" --snplist "$COVERAGE_SNPLIST" \
    --out "$CLEAN_SNPLIST"

  echo "=== Track $TRACK: shared-matrix fixed reference on the cleaned marker set (09) ==="
  python3 scripts/ecotype_pca_v2/09_export_fixed_reference_eigenstrat.py \
    --panel C --library-type pooled_mixed --track "$TRACK" \
    --panel-snp "$PANEL_SNP" --panel-geno "$PANEL_GENO" --panel-ind "$PANEL_IND" \
    --fixed-snplist "$CLEAN_SNPLIST" --reference-keep "$CIVAN_REFERENCE_KEEP" \
    --label civan --out-dir "$TRACK_DIR"
  REFERENCE_PREFIX="$TRACK_DIR/civan.pooled_mixed.$TRACK.fixed_reference"

  CALLS_ARGS=()
  for SAMPLE in $CIVAN_ANCIENT_SAMPLES; do
    BAM="$BAMDIR/$SAMPLE.besthit_oryza.irgsp.bam"
    if [[ ! -s "$BAM" ]]; then
      echo "WARNING: $TRACK/$SAMPLE: BAM missing or empty, technical_execution=FAIL" >&2
      printf 'shared\t%s\t%s\tFAIL\tBAM missing or empty\tNA\n' "$TRACK" "$SAMPLE" >> "$SUMMARY"
      continue
    fi
    set +e
    python3 scripts/ecotype_pca_v2/10_call_ancient_fixed_markers.py \
      --config "$RICE_PCA_CONFIG" --bam "$BAM" --fixed-snp "$REFERENCE_PREFIX.snp" \
      --sample "$SAMPLE" --panel C --library-type pooled_mixed --track "$TRACK" \
      --out-dir "$TRACK_DIR/calls" 2> "$TRACK_DIR/calls/$SAMPLE.$TRACK.stderr.log"
    RC=$?
    set -e
    if [[ $RC -eq 0 ]]; then
      printf 'shared\t%s\t%s\tPASS\t\tNA\n' "$TRACK" "$SAMPLE" >> "$SUMMARY"
      CALLS_ARGS+=("$SAMPLE=$TRACK_DIR/calls/$SAMPLE.C.pooled_mixed.$TRACK.calls.txt")
    else
      NOTE=$(tail -1 "$TRACK_DIR/calls/$SAMPLE.$TRACK.stderr.log" 2>/dev/null || echo "unknown error")
      echo "WARNING: $TRACK/$SAMPLE: calling failed (exit $RC): $NOTE" >&2
      printf 'shared\t%s\t%s\tFAIL\t%s\tNA\n' "$TRACK" "$SAMPLE" "$NOTE" >> "$SUMMARY"
    fi
  done
  if [[ ${#CALLS_ARGS[@]} -eq 0 ]]; then
    echo "FATAL: track $TRACK: no sample completed calling against the MAF/LD-cleaned marker set" >&2
    exit 3
  fi

  python3 scripts/ecotype_pca_v2/11_build_ancient_callability.py \
    --config "$RICE_PCA_CONFIG" --fixed-snp "$REFERENCE_PREFIX.snp" \
    $(for a in "${CALLS_ARGS[@]}"; do printf -- '--calls %s ' "$a"; done) \
    --panel C --library-type pooled_mixed --track "$TRACK" \
    --out "$TRACK_DIR/civan.$TRACK.callability.tsv"

  python3 scripts/ecotype_pca_v2/13_merge_ancients_fixed_panel.py \
    --reference-geno "$REFERENCE_PREFIX.eigenstratgeno" --reference-ind "$REFERENCE_PREFIX.ind" \
    --fixed-snp "$REFERENCE_PREFIX.snp" \
    $(for a in "${CALLS_ARGS[@]}"; do printf -- '--calls %s ' "$a"; done) \
    --ancient-poplabel Ancient --label civan --out-dir "$TRACK_DIR"

  bash scripts/ecotype_pca_v2/14_run_fixed_smartpca.sh \
    --config "$RICE_PCA_CONFIG" \
    --geno "$TRACK_DIR/civan.merged.eigenstratgeno" --snp "$REFERENCE_PREFIX.snp" \
    --ind "$TRACK_DIR/civan.merged.ind" --poplist "$REFERENCE_PREFIX.poplistname" \
    --label "civan.$TRACK.pca" --out-dir "$TRACK_DIR"

  EXPECTED_N=$(wc -l < "$TRACK_DIR/civan.merged.ind")
  python3 scripts/ecotype_pca_v2/15_pca_qc.py \
    --evec "$TRACK_DIR/civan.$TRACK.pca.evec" --ind "$TRACK_DIR/civan.merged.ind" \
    --expected-n "$EXPECTED_N" --out "$TRACK_DIR/civan.$TRACK.pca_qc.tsv"

  python3 scripts/ecotype_pca_v2/16_projection_summary.py \
    --panel C --evec "$TRACK_DIR/civan.$TRACK.pca.evec" \
    --out "$TRACK_DIR/civan.$TRACK.projection_summary.tsv"

  REPORT_ARGS=()
  for a in "${CALLS_ARGS[@]}"; do
    SAMPLE="${a%%=*}"
    REPORT_ARGS+=("$SAMPLE=$TRACK_DIR/calls/$SAMPLE.C.pooled_mixed.$TRACK.call_report.tsv")
  done
  python3 scripts/ecotype_pca_v2/22_classify_scientific_projection.py \
    $(for a in "${REPORT_ARGS[@]}"; do printf -- '--call-report %s ' "$a"; done) \
    --out "$TRACK_DIR/civan.$TRACK.scientific_projection.tsv"

  bash scripts/ecotype_pca_v2/26_plot_pc_pairs.sh \
    --evec-label "civan_${TRACK}_shared" --evec "$TRACK_DIR/civan.$TRACK.pca.evec" \
    --title-prefix "Civan MAF/LD-cleaned shared matrix ($TRACK track)" \
    --out-prefix "$TRACK_DIR/civan.$TRACK.shared"

  echo "=== Track $TRACK: per-sample private-axis projections (v1-style, own covered sites only) ==="
  PRIVATE_TRACK_DIR="$OUT/PRIVATE/$TRACK"
  mkdir -p "$PRIVATE_TRACK_DIR"
  for a in "${CALLS_ARGS[@]}"; do
    SAMPLE="${a%%=*}"
    SDIR="$PRIVATE_TRACK_DIR/$SAMPLE"
    mkdir -p "$SDIR/calls"
    CALL_SITES="$TRACK_DIR/calls/$SAMPLE.C.pooled_mixed.$TRACK.call_sites.tsv"
    PRIVATE_SNPLIST="$SDIR/${SAMPLE}.private.snplist"

    set +e
    python3 scripts/ecotype_pca_v2/24_extract_sample_covered_sites.py \
      --call-sites "$CALL_SITES" --out "$PRIVATE_SNPLIST"
    RC=$?
    set -e
    if [[ $RC -ne 0 ]]; then
      echo "WARNING: $TRACK/$SAMPLE: zero covered sites in the MAF/LD-cleaned set, skipping private axis" >&2
      printf 'private\t%s\t%s\tFAIL\tzero covered sites in cleaned marker set\t0\n' "$TRACK" "$SAMPLE" >> "$SUMMARY"
      continue
    fi

    python3 scripts/ecotype_pca_v2/09_export_fixed_reference_eigenstrat.py \
      --panel C --library-type pooled_mixed --track "$TRACK" \
      --panel-snp "$PANEL_SNP" --panel-geno "$PANEL_GENO" --panel-ind "$PANEL_IND" \
      --fixed-snplist "$PRIVATE_SNPLIST" --reference-keep "$CIVAN_REFERENCE_KEEP" \
      --label "civan.${SAMPLE}" --out-dir "$SDIR"
    PRIVATE_REF_PREFIX="$SDIR/civan.${SAMPLE}.pooled_mixed.$TRACK.fixed_reference"

    set +e
    python3 scripts/ecotype_pca_v2/10_call_ancient_fixed_markers.py \
      --config "$RICE_PCA_CONFIG" --bam "$BAMDIR/$SAMPLE.besthit_oryza.irgsp.bam" \
      --fixed-snp "$PRIVATE_REF_PREFIX.snp" --sample "$SAMPLE" --panel C \
      --library-type pooled_mixed --track "$TRACK" --out-dir "$SDIR/calls" \
      2> "$SDIR/calls/$SAMPLE.$TRACK.private.stderr.log"
    RC=$?
    set -e
    if [[ $RC -ne 0 ]]; then
      NOTE=$(tail -1 "$SDIR/calls/$SAMPLE.$TRACK.private.stderr.log" 2>/dev/null || echo "unknown error")
      echo "WARNING: $TRACK/$SAMPLE: private-axis calling failed: $NOTE" >&2
      printf 'private\t%s\t%s\tFAIL\t%s\tNA\n' "$TRACK" "$SAMPLE" "$NOTE" >> "$SUMMARY"
      continue
    fi
    PRIVATE_CALLABLE_N=$(grep -cE '^[02]$' "$SDIR/calls/$SAMPLE.C.pooled_mixed.$TRACK.calls.txt" || true)

    python3 scripts/ecotype_pca_v2/13_merge_ancients_fixed_panel.py \
      --reference-geno "$PRIVATE_REF_PREFIX.eigenstratgeno" --reference-ind "$PRIVATE_REF_PREFIX.ind" \
      --fixed-snp "$PRIVATE_REF_PREFIX.snp" \
      --calls "$SAMPLE=$SDIR/calls/$SAMPLE.C.pooled_mixed.$TRACK.calls.txt" \
      --ancient-poplabel Ancient --label "civan.${SAMPLE}" --out-dir "$SDIR"

    bash scripts/ecotype_pca_v2/14_run_fixed_smartpca.sh \
      --config "$RICE_PCA_CONFIG" \
      --geno "$SDIR/civan.${SAMPLE}.merged.eigenstratgeno" --snp "$PRIVATE_REF_PREFIX.snp" \
      --ind "$SDIR/civan.${SAMPLE}.merged.ind" --poplist "$PRIVATE_REF_PREFIX.poplistname" \
      --label "civan.${SAMPLE}.$TRACK.pca" --out-dir "$SDIR"

    PRIVATE_EXPECTED_N=$(wc -l < "$SDIR/civan.${SAMPLE}.merged.ind")
    python3 scripts/ecotype_pca_v2/15_pca_qc.py \
      --evec "$SDIR/civan.${SAMPLE}.$TRACK.pca.evec" --ind "$SDIR/civan.${SAMPLE}.merged.ind" \
      --expected-n "$PRIVATE_EXPECTED_N" --out "$SDIR/civan.${SAMPLE}.$TRACK.pca_qc.tsv"

    python3 scripts/ecotype_pca_v2/16_projection_summary.py \
      --panel C --evec "$SDIR/civan.${SAMPLE}.$TRACK.pca.evec" \
      --out "$SDIR/civan.${SAMPLE}.$TRACK.projection_summary.tsv"

    bash scripts/ecotype_pca_v2/26_plot_pc_pairs.sh \
      --evec-label "civan_${TRACK}_${SAMPLE}_private" --evec "$SDIR/civan.${SAMPLE}.$TRACK.pca.evec" \
      --title-prefix "Civan private axis: $SAMPLE ($TRACK track, own covered sites only)" \
      --out-prefix "$SDIR/civan.${SAMPLE}.$TRACK.private"

    printf 'private\t%s\t%s\tPASS\t\t%s\n' "$TRACK" "$SAMPLE" "$PRIVATE_CALLABLE_N" >> "$SUMMARY"
    echo "PASS: $TRACK/$SAMPLE private axis complete (callable_n=$PRIVATE_CALLABLE_N)"
  done
done

echo "PASS: 51_civan_maf_ld_and_private_axis complete"
echo "See $SUMMARY, SHARED/*/civan.*.scientific_projection.tsv, and */*.png"
