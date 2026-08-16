#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?stage 50 must run in SLURM}"
: "${RICE_PCA_REPO_ROOT:?}"; : "${RICE_PCA_CONFIG:?}"; : "${RICE_PCA_ATTEMPT_DIR:?}"
: "${CIVAN_UNION_SITES:?path to ancient_union_sites.tsv from 19_survey_ancient_coverage.py}"
: "${CIVAN_UNION_SITES_TV:?path to ancient_union_sites.TV.tsv from 20_filter_coverage_sites_to_transversions.py}"
: "${CIVAN_REFERENCE_KEEP:?panel C reference/axis-builder keep-list -- the output of \
06_build_reference_sample_set.py --panel C (a *.reference_samples.keep file). Locate an \
existing one with: find \$RESULTS_V2_ROOT -name '*.reference_samples.keep'. If none exists \
yet, run 06_build_reference_sample_set.py --panel C --ind-file <civan .filtered.ind> first \
-- do not guess this path.}"
: "${CIVAN_ANCIENT_SAMPLES:?space-separated list of ancient sample IDs to project, e.g. \
all 16 besthit samples}"
cd "$RICE_PCA_REPO_ROOT"

eval "$(python3 - "$RICE_PCA_CONFIG" <<'PY'
import shlex,sys,yaml
c=yaml.safe_load(open(sys.argv[1])); p=c['inputs']['panel_C_civan']
for k,v in {
  'PANEL_SNP': f"{p['dir']}/{p['prefix']}.snp",
  'PANEL_GENO': f"{p['dir']}/{p['prefix']}{p['filtered_suffix']}.eigenstratgeno",
  'PANEL_IND': f"{p['dir']}/{p['prefix']}{p['filtered_suffix']}.ind",
  'BAMDIR': c['inputs']['ancient_bam_dir'],
}.items(): print(f'{k}={shlex.quote(v)}')
PY
)"

OUT="$RICE_PCA_ATTEMPT_DIR"
SUMMARY="$OUT/stage50_summary.tsv"
printf 'track\tsample\tbam_path\ttechnical_execution\ttechnical_note\n' > "$SUMMARY"

read -r -a ANCIENT_SAMPLE_IDS <<< "$CIVAN_ANCIENT_SAMPLES"
[[ ${#ANCIENT_SAMPLE_IDS[@]} -gt 0 ]] || { echo "FATAL: no ancient samples supplied" >&2; exit 2; }
SEEN_SAMPLE_IDS=" "
for SAMPLE in "${ANCIENT_SAMPLE_IDS[@]}"; do
  [[ "$SAMPLE" =~ ^[A-Za-z0-9._-]+$ ]] || { echo "FATAL: invalid ancient sample ID: $SAMPLE" >&2; exit 2; }
  [[ "$SEEN_SAMPLE_IDS" != *" $SAMPLE "* ]] || { echo "FATAL: duplicate ancient sample ID: $SAMPLE" >&2; exit 2; }
  SEEN_SAMPLE_IDS+="$SAMPLE "
done

for TRACK in ALL TV; do
  echo "=== track $TRACK ==="
  case "$TRACK" in
    ALL) SITES="$CIVAN_UNION_SITES" ;;
    TV)  SITES="$CIVAN_UNION_SITES_TV" ;;
  esac
  TRACK_DIR="$OUT/$TRACK"
  mkdir -p "$TRACK_DIR/calls"

  FIXED_LIST="$TRACK_DIR/fixed_snplist.txt"
  python3 scripts/ecotype_pca_v2/21_extract_fixed_snplist.py \
    --sites-tsv "$SITES" --out "$FIXED_LIST"

  python3 scripts/ecotype_pca_v2/09_export_fixed_reference_eigenstrat.py \
    --panel C --library-type pooled_mixed --track "$TRACK" \
    --panel-snp "$PANEL_SNP" --panel-geno "$PANEL_GENO" --panel-ind "$PANEL_IND" \
    --fixed-snplist "$FIXED_LIST" --reference-keep "$CIVAN_REFERENCE_KEEP" \
    --label civan --out-dir "$TRACK_DIR"

  REFERENCE_PREFIX="$TRACK_DIR/civan.pooled_mixed.$TRACK.fixed_reference"

  CALLS_ARGS=()
  REPORT_ARGS=()
  TECHNICAL_FAILURE_N=0
  for SAMPLE in "${ANCIENT_SAMPLE_IDS[@]}"; do
    BAM="$BAMDIR/$SAMPLE.besthit_oryza.irgsp.bam"
    REPORT_ARGS+=("$SAMPLE=$TRACK_DIR/calls/$SAMPLE.C.pooled_mixed.$TRACK.call_report.tsv")
    if [[ ! -s "$BAM" ]]; then
      echo "WARNING: $TRACK/$SAMPLE: BAM missing or empty, technical_execution=FAIL" >&2
      printf '%s\t%s\t%s\tFAIL\t%s\n' "$TRACK" "$SAMPLE" "$BAM" "BAM missing or empty" >> "$SUMMARY"
      TECHNICAL_FAILURE_N=$((TECHNICAL_FAILURE_N + 1))
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
      printf '%s\t%s\t%s\tPASS\t\n' "$TRACK" "$SAMPLE" "$BAM" >> "$SUMMARY"
      CALLS_ARGS+=("$SAMPLE=$TRACK_DIR/calls/$SAMPLE.C.pooled_mixed.$TRACK.calls.txt")
    else
      NOTE=$(tail -1 "$TRACK_DIR/calls/$SAMPLE.$TRACK.stderr.log" 2>/dev/null || echo "unknown error")
      echo "WARNING: $TRACK/$SAMPLE: calling step failed (exit $RC): $NOTE" >&2
      printf '%s\t%s\t%s\tFAIL\t%s\n' "$TRACK" "$SAMPLE" "$BAM" "$NOTE" >> "$SUMMARY"
      TECHNICAL_FAILURE_N=$((TECHNICAL_FAILURE_N + 1))
    fi
  done

  REPORT_CLI=()
  for a in "${REPORT_ARGS[@]}"; do REPORT_CLI+=(--call-report "$a"); done
  python3 scripts/ecotype_pca_v2/22_classify_scientific_projection.py \
    "${REPORT_CLI[@]}" \
    --out "$TRACK_DIR/civan.$TRACK.scientific_projection.tsv"

  if [[ ${#CALLS_ARGS[@]} -eq 0 ]]; then
    echo "FATAL: track $TRACK: no sample completed the calling step" >&2
    exit 3
  fi
  if [[ $TECHNICAL_FAILURE_N -ne 0 ]]; then
    echo "FATAL: track $TRACK: $TECHNICAL_FAILURE_N/${#ANCIENT_SAMPLE_IDS[@]} sample(s) failed technical calling; refusing an incomplete Stage 50 receipt" >&2
    exit 3
  fi

  CALLS_CLI=()
  for a in "${CALLS_ARGS[@]}"; do CALLS_CLI+=(--calls "$a"); done
  python3 scripts/ecotype_pca_v2/11_build_ancient_callability.py \
    --config "$RICE_PCA_CONFIG" --fixed-snp "$REFERENCE_PREFIX.snp" \
    "${CALLS_CLI[@]}" \
    --panel C --library-type pooled_mixed --track "$TRACK" \
    --out "$TRACK_DIR/civan.$TRACK.callability.tsv"

  python3 scripts/ecotype_pca_v2/13_merge_ancients_fixed_panel.py \
    --reference-geno "$REFERENCE_PREFIX.eigenstratgeno" --reference-ind "$REFERENCE_PREFIX.ind" \
    --fixed-snp "$REFERENCE_PREFIX.snp" \
    "${CALLS_CLI[@]}" \
    --ancient-poplabel Ancient --label civan --out-dir "$TRACK_DIR"

  bash scripts/ecotype_pca_v2/14_run_fixed_smartpca.sh \
    --config "$RICE_PCA_CONFIG" \
    --geno "$TRACK_DIR/civan.merged.eigenstratgeno" \
    --snp "$REFERENCE_PREFIX.snp" \
    --ind "$TRACK_DIR/civan.merged.ind" \
    --poplist "$REFERENCE_PREFIX.poplistname" \
    --label "civan.$TRACK.pca" --out-dir "$TRACK_DIR"

  EXPECTED_N=$(( $(wc -l < "$TRACK_DIR/civan.merged.ind") ))
  python3 scripts/ecotype_pca_v2/15_pca_qc.py \
    --evec "$TRACK_DIR/civan.$TRACK.pca.evec" --ind "$TRACK_DIR/civan.merged.ind" \
    --expected-n "$EXPECTED_N" --projection-label Ancient \
    --out "$TRACK_DIR/civan.$TRACK.pca_qc.tsv"

  python3 scripts/ecotype_pca_v2/16_projection_summary.py \
    --panel C --evec "$TRACK_DIR/civan.$TRACK.pca.evec" \
    --out "$TRACK_DIR/civan.$TRACK.projection_summary.tsv"

  echo "PASS: track $TRACK complete ($(( ${#CALLS_ARGS[@]} )) sample(s) technically succeeded)"
done

echo "PASS: 50_civan_coverage_aware_projection -- technical_execution and scientific_projection recorded per sample per track in */civan.*.scientific_projection.tsv and $SUMMARY"
echo "NOTE: this stage does not auto-unlock Stage 60 -- unlock conditions are unchanged"
