#!/usr/bin/env bash
# test_integration_synthetic.sh
#
# End-to-end integration tests against a tiny SYNTHETIC panel (not real
# server data) exercising plink2/convertf/bedtools the way 00/01/02/06/07/08
# actually call them. Item 15 of the Batch 1 correction request.
#
# NOT executed by Claude Code -- no server access (github-repo-protocol
# Rule 3). Written and reviewed carefully, but the pass/fail transcript in
# the correction-batch report is from test_lib_ecotype_v2.py,
# test_08_streaming.py and test_06_and_01_logic.py only (pure-python, no
# external tools, actually run locally). This script needs plink2,
# convertf, and bedtools on PATH -- run it on the server and paste back the
# output; that's the only way its results become real rather than assumed.
#
# Usage: bash test_integration_synthetic.sh /path/to/scripts/ecotype_pca_v2
#
# Builds everything under a fresh temp dir and cleans up on exit (successful
# or not) unless KEEP_TEST_DIR=1 is set in the environment.

set -uo pipefail  # NOT -e: we want to keep running and tally pass/fail

SCRIPTS_DIR="${1:?Usage: $0 /path/to/scripts/ecotype_pca_v2}"
CONFIG="$SCRIPTS_DIR/config/ecotype_pca_v2.yaml"
[[ -f "$CONFIG" ]] || { echo "FATAL: config not found at $CONFIG"; exit 1; }

TMP="$(mktemp -d)"
[[ "${KEEP_TEST_DIR:-0}" == "1" ]] || trap 'rm -rf "$TMP"' EXIT
echo "test workspace: $TMP"

PASS=0
FAIL=0
FAILED_NAMES=()

expect_exit() {
  local name="$1" expected="$2"; shift 2
  set +e
  "$@" > "$TMP/${name}.out" 2>&1
  local actual=$?
  set -e
  if [[ "$actual" == "$expected" ]]; then
    echo "PASS: $name (exit $actual, expected $expected)"
    PASS=$((PASS+1))
  else
    echo "FAIL: $name (exit $actual, expected $expected) -- see $TMP/${name}.out"
    tail -20 "$TMP/${name}.out" | sed 's/^/    /'
    FAIL=$((FAIL+1))
    FAILED_NAMES+=("$name")
  fi
}

expect_file_exists() {
  local name="$1" path="$2"
  if [[ -f "$path" ]]; then
    echo "PASS: $name (file exists: $path)"
    PASS=$((PASS+1))
  else
    echo "FAIL: $name (file MISSING: $path)"
    FAIL=$((FAIL+1))
    FAILED_NAMES+=("$name")
  fi
}

expect_line_count() {
  local name="$1" path="$2" expected="$3"
  local actual
  actual=$(wc -l < "$path" 2>/dev/null || echo -1)
  if [[ "$actual" == "$expected" ]]; then
    echo "PASS: $name ($path has $actual lines, expected $expected)"
    PASS=$((PASS+1))
  else
    echo "FAIL: $name ($path has $actual lines, expected $expected)"
    FAIL=$((FAIL+1))
    FAILED_NAMES+=("$name")
  fi
}

# ============================================================
# 1. --help exits 0 on every script (item 13)
# ============================================================
for s in 00_validate_inputs.py 01_make_panel_manifest.py 03_audit_panel.py \
         04_audit_720_ld.py 05_intersect_panel_baits.py 06_build_reference_sample_set.py \
         08_make_5kb_thinned_markers.py; do
  expect_exit "help_exit0_${s}" 0 python3 "$SCRIPTS_DIR/$s" --help
done
for s in 02_convert_eigenstrat_for_plink.sh 07_make_fixed_markers.sh; do
  expect_exit "help_exit0_${s}" 0 bash "$SCRIPTS_DIR/$s" --help
done

# ============================================================
# 2. build a tiny synthetic Panel-A-like EIGENSTRAT panel
#    12 samples (IND/AUS/ARO/TRJ/TEJ x a few each), 20 SNPs on chr1/chr2,
#    mix of TV and transition pairs, all biallelic ACGT, unique IDs/pos.
# ============================================================
PANEL_DIR="$TMP/panel_synth"
mkdir -p "$PANEL_DIR"
PREFIX="synthA"

cat > "$PANEL_DIR/${PREFIX}.ind" <<EOF
s01	U	IND
s02	U	IND
s03	U	AUS
s04	U	AUS
s05	U	ARO
s06	U	ARO
s07	U	TRJ
s08	U	TRJ
s09	U	TEJ
s10	U	TEJ
s11	U	ADM
s12	U	Ancient
EOF
cp "$PANEL_DIR/${PREFIX}.ind" "$PANEL_DIR/${PREFIX}.filtered.ind"

# 20 SNPs: alternate TV (A/C) and transition (A/G) alleles, unique pos, chr1+chr2
python3 - "$PANEL_DIR/${PREFIX}.snp" <<'PYEOF'
import sys
path = sys.argv[1]
tv_pairs = [("A","C"), ("A","T"), ("C","G"), ("G","T")]
ts_pairs = [("A","G"), ("C","T")]
with open(path, "w") as fh:
    pos = 1000
    for i in range(20):
        chrom = "1" if i < 12 else "2"
        pair = tv_pairs[i % 4] if i % 2 == 0 else ts_pairs[i % 2]
        fh.write(f"snp{i+1}\t{chrom}\t0\t{pos}\t{pair[0]}\t{pair[1]}\n")
        pos += 1000
PYEOF

# genotypes: 20 SNPs x 12 samples, EIGENSTRAT text (rows=SNPs, cols=samples),
# random-ish but deterministic 0/1/2, no missing (9) so geno/MAF filters have
# something nonzero to work with.
python3 - "$PANEL_DIR/${PREFIX}.eigenstratgeno" <<'PYEOF'
import sys, random
path = sys.argv[1]
rng = random.Random(1)
with open(path, "w") as fh:
    for _ in range(20):
        fh.write("".join(str(rng.choice([0,1,2])) for _ in range(12)) + "\n")
PYEOF
cp "$PANEL_DIR/${PREFIX}.eigenstratgeno" "$PANEL_DIR/${PREFIX}.filtered.eigenstratgeno"

# ============================================================
# 3. 02: EIGENSTRAT -> PLINK
# ============================================================
CONVERT_OUT="$TMP/convert"
expect_exit "02_convert_synthA" 0 bash "$SCRIPTS_DIR/02_convert_eigenstrat_for_plink.sh" \
  --dir "$PANEL_DIR" --prefix "$PREFIX" --out-dir "$CONVERT_OUT"
BFILE="$CONVERT_OUT/${PREFIX}.plink"
expect_file_exists "02_output_bim" "${BFILE}.bim"
expect_file_exists "02_output_fam" "${BFILE}.fam"

# ============================================================
# 4. 06: reference sample set, all 5 labels present -> should succeed
# ============================================================
REF_OUT="$TMP/ref"
expect_exit "06_panel_A_all_labels_present" 0 python3 "$SCRIPTS_DIR/06_build_reference_sample_set.py" \
  --config "$CONFIG" --panel A --label synthA \
  --ind-file "$PANEL_DIR/${PREFIX}.filtered.ind" --fam-file "${BFILE}.fam" \
  --out-dir "$REF_OUT"
expect_line_count "06_keep_list_has_10_reference_samples" "$REF_OUT/synthA.reference_samples.keep" 10

# 06: missing label (drop TEJ from a copy of the .ind) -> must hard-fail (item 11)
MISSING_LABEL_IND="$TMP/synthA_missing_tej.ind"
grep -v "TEJ" "$PANEL_DIR/${PREFIX}.filtered.ind" > "$MISSING_LABEL_IND"
expect_exit "06_panel_A_missing_label_hard_fails" 3 python3 "$SCRIPTS_DIR/06_build_reference_sample_set.py" \
  --config "$CONFIG" --panel A --label synthA_missing \
  --ind-file "$MISSING_LABEL_IND" --fam-file "${BFILE}.fam" \
  --out-dir "$TMP/ref_missing"

# 06: duplicate sample ID -> must hard-fail (item 11)
DUP_IND="$TMP/synthA_dup.ind"
{ cat "$PANEL_DIR/${PREFIX}.filtered.ind"; echo -e "s01\tU\tAUS"; } > "$DUP_IND"
expect_exit "06_duplicate_sample_id_hard_fails" 3 python3 "$SCRIPTS_DIR/06_build_reference_sample_set.py" \
  --config "$CONFIG" --panel A --label synthA_dup \
  --ind-file "$DUP_IND" --fam-file "${BFILE}.fam" \
  --out-dir "$TMP/ref_dup"

# ============================================================
# 5. 07: enum guards (item 7) -- every one of these must hard-fail, exit 1
# ============================================================
FIXED_OUT="$TMP/fixed"
expect_exit "07_bad_panel_enum" 1 bash "$SCRIPTS_DIR/07_make_fixed_markers.sh" \
  --config "$CONFIG" --panel Z --sensitivity primary --library-type shotgun --track TV \
  --bfile "$BFILE" --keep "$REF_OUT/synthA.reference_samples.keep" --label synthA --out-dir "$FIXED_OUT"
expect_exit "07_bad_library_type_enum" 1 bash "$SCRIPTS_DIR/07_make_fixed_markers.sh" \
  --config "$CONFIG" --panel A --sensitivity primary --library-type wgs --track TV \
  --bfile "$BFILE" --keep "$REF_OUT/synthA.reference_samples.keep" --label synthA --out-dir "$FIXED_OUT"
expect_exit "07_bad_track_enum" 1 bash "$SCRIPTS_DIR/07_make_fixed_markers.sh" \
  --config "$CONFIG" --panel A --sensitivity primary --library-type shotgun --track transversions \
  --bfile "$BFILE" --keep "$REF_OUT/synthA.reference_samples.keep" --label synthA --out-dir "$FIXED_OUT"
expect_exit "07_bad_sensitivity_enum" 1 bash "$SCRIPTS_DIR/07_make_fixed_markers.sh" \
  --config "$CONFIG" --panel A --sensitivity thinning_only --library-type shotgun --track TV \
  --bfile "$BFILE" --keep "$REF_OUT/synthA.reference_samples.keep" --label synthA --out-dir "$FIXED_OUT"
expect_exit "07_capture_without_snp_list_fails" 1 bash "$SCRIPTS_DIR/07_make_fixed_markers.sh" \
  --config "$CONFIG" --panel A --sensitivity primary --library-type capture --track TV \
  --bfile "$BFILE" --keep "$REF_OUT/synthA.reference_samples.keep" --label synthA --out-dir "$FIXED_OUT"

# ============================================================
# 6. 07: real runs -- TV vs ALL, primary vs S1-S4, geno_maf_only stage
# ============================================================
expect_exit "07_shotgun_TV_primary" 0 bash "$SCRIPTS_DIR/07_make_fixed_markers.sh" \
  --config "$CONFIG" --panel A --sensitivity primary --library-type shotgun --track TV \
  --bfile "$BFILE" --keep "$REF_OUT/synthA.reference_samples.keep" --label synthA --out-dir "$FIXED_OUT"
expect_file_exists "07_TV_primary_fixed_snplist" "$FIXED_OUT/synthA.shotgun.TV.primary.fixed.snplist"
expect_file_exists "07_TV_primary_manifest_has_split_fields" "$FIXED_OUT/synthA.shotgun.TV.primary.marker_manifest.tsv"
if [[ -f "$FIXED_OUT/synthA.shotgun.TV.primary.marker_manifest.tsv" ]]; then
  if head -1 "$FIXED_OUT/synthA.shotgun.TV.primary.marker_manifest.tsv" | grep -q "after_site_missingness" && \
     head -1 "$FIXED_OUT/synthA.shotgun.TV.primary.marker_manifest.tsv" | grep -q "after_MAF"; then
    echo "PASS: 07_manifest_has_separate_site_missingness_and_MAF_columns (item 9)"
    PASS=$((PASS+1))
  else
    echo "FAIL: 07_manifest_has_separate_site_missingness_and_MAF_columns (item 9)"
    FAIL=$((FAIL+1)); FAILED_NAMES+=("07_manifest_split_fields")
  fi
fi

expect_exit "07_shotgun_ALL_S1" 0 bash "$SCRIPTS_DIR/07_make_fixed_markers.sh" \
  --config "$CONFIG" --panel A --sensitivity S1 --library-type shotgun --track ALL \
  --bfile "$BFILE" --keep "$REF_OUT/synthA.reference_samples.keep" --label synthA --out-dir "$FIXED_OUT"
expect_exit "07_shotgun_TV_S4" 0 bash "$SCRIPTS_DIR/07_make_fixed_markers.sh" \
  --config "$CONFIG" --panel A --sensitivity S4 --library-type shotgun --track TV \
  --bfile "$BFILE" --keep "$REF_OUT/synthA.reference_samples.keep" --label synthA --out-dir "$FIXED_OUT"

# geno_maf_only stage (item 1: this is the interface that was broken before)
GENOMAF_OUT="$TMP/genomaf"
expect_exit "07_geno_maf_only_stage" 0 bash "$SCRIPTS_DIR/07_make_fixed_markers.sh" \
  --config "$CONFIG" --panel A --sensitivity primary --library-type shotgun --track ALL \
  --stage geno_maf_only \
  --bfile "$BFILE" --keep "$REF_OUT/synthA.reference_samples.keep" --label synthA --out-dir "$GENOMAF_OUT"
expect_file_exists "07_geno_maf_only_bim" "$GENOMAF_OUT/synthA.shotgun.ALL.primary.geno_maf_filtered.bim"
expect_file_exists "07_geno_maf_only_manifest" "$GENOMAF_OUT/synthA.shotgun.ALL.primary.geno_maf_manifest.tsv"
# and the LD-pruned .fixed.snplist must NOT exist for this run (stage stopped early)
if [[ ! -e "$GENOMAF_OUT/synthA.shotgun.ALL.primary.fixed.snplist" ]]; then
  echo "PASS: 07_geno_maf_only_does_not_produce_fixed_snplist"
  PASS=$((PASS+1))
else
  echo "FAIL: 07_geno_maf_only_does_not_produce_fixed_snplist (should not exist yet)"
  FAIL=$((FAIL+1)); FAILED_NAMES+=("07_geno_maf_only_no_fixed")
fi

# overwrite protection (item 14): re-running without --overwrite must fail
expect_exit "07_overwrite_protection" 3 bash "$SCRIPTS_DIR/07_make_fixed_markers.sh" \
  --config "$CONFIG" --panel A --sensitivity primary --library-type shotgun --track TV \
  --bfile "$BFILE" --keep "$REF_OUT/synthA.reference_samples.keep" --label synthA --out-dir "$FIXED_OUT"

# ============================================================
# 7. 08: paperlike_5kb thinning against 07's geno_maf_only output
# ============================================================
THIN_OUT="$TMP/thin"
expect_exit "08_paperlike_5kb" 0 python3 "$SCRIPTS_DIR/08_make_5kb_thinned_markers.py" \
  --config "$CONFIG" --label synthA \
  --geno-maf-bim "$GENOMAF_OUT/synthA.shotgun.ALL.primary.geno_maf_filtered.bim" \
  --upstream-manifest "$GENOMAF_OUT/synthA.shotgun.ALL.primary.geno_maf_manifest.tsv" \
  --out-dir "$THIN_OUT"
expect_file_exists "08_output_snplist" "$THIN_OUT/synthA.paperlike_5kb.fixed.snplist"
expect_file_exists "08_output_manifest" "$THIN_OUT/synthA.paperlike_5kb.marker_manifest.tsv"

# ============================================================
# 8. summary
# ============================================================
echo
echo "============================================================"
echo "TOTAL: $PASS passed, $FAIL failed"
if [[ $FAIL -gt 0 ]]; then
  echo "FAILED TESTS: ${FAILED_NAMES[*]}"
  exit 1
fi
exit 0
