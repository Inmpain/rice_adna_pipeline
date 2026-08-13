#!/usr/bin/env bash
# Phase 1 / Batch 1, script 07.
#
# Freezes ONE coordinate system's marker set: reference samples only ->
# (capture bait intersect, if --library-type capture) -> TV/ALL split ->
# site missingness -> MAF -> LD prune (or skip, for callers that only want
# the pre-LD intermediate, e.g. 08's paperlike_5kb route) -> *.fixed.snplist
# + *.marker_manifest.tsv + md5.
#
# All numeric parameters (geno/maf/ld window/r2) are resolved from config.yaml
# via lib_ecotype_v2.resolve_marker_params -- this script never accepts a raw
# --maf / --geno / --ld-window-kb override on the command line, by design,
# so a debugging session cannot accidentally drift a frozen parameter.
#
# Once *.fixed.snplist exists for a given (panel, library_type, track,
# sensitivity) combination, it must never be regenerated to add/drop SNPs
# because an ancient sample was added -- pass --overwrite only when you mean
# to intentionally rebuild (e.g. after fixing an implementation bug in this
# script itself, not to chase a different SNP count).

set -euo pipefail

usage() {
  cat <<EOF
Usage: $0 --config CFG --panel {A|B|C} --sensitivity {primary|S1|S2|S3|S4} \\
          --library-type {shotgun|capture} --track {TV|ALL} \\
          --bfile PLINK_PREFIX --keep REFERENCE_KEEP_FILE --label LABEL \\
          --out-dir OUT_DIR [--capture-snp-list FILE] [--overwrite]

  --capture-snp-list   required if --library-type capture: output of
                        05_intersect_panel_baits.py ({label}.capture_compatible.snp)
EOF
  exit 1
}

CONFIG=""; PANEL=""; SENS=""; LIBTYPE=""; TRACK=""; BFILE=""; KEEP=""; LABEL=""; OUT_DIR=""
CAPTURE_SNP_LIST=""; OVERWRITE=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --config) CONFIG="$2"; shift 2 ;;
    --panel) PANEL="$2"; shift 2 ;;
    --sensitivity) SENS="$2"; shift 2 ;;
    --library-type) LIBTYPE="$2"; shift 2 ;;
    --track) TRACK="$2"; shift 2 ;;
    --bfile) BFILE="$2"; shift 2 ;;
    --keep) KEEP="$2"; shift 2 ;;
    --label) LABEL="$2"; shift 2 ;;
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    --capture-snp-list) CAPTURE_SNP_LIST="$2"; shift 2 ;;
    --overwrite) OVERWRITE=1; shift ;;
    -h|--help) usage ;;
    *) echo "unknown argument: $1" >&2; usage ;;
  esac
done
for v in CONFIG PANEL SENS LIBTYPE TRACK BFILE KEEP LABEL OUT_DIR; do
  [[ -z "${!v}" ]] && { echo "missing required --${v,,}" >&2; usage; }
done
if [[ "$LIBTYPE" == "capture" && -z "$CAPTURE_SNP_LIST" ]]; then
  echo "FATAL: --library-type capture requires --capture-snp-list (05's output)" >&2
  exit 1
fi

mkdir -p "$OUT_DIR"
RUN_ID="${LABEL}.${LIBTYPE}.${TRACK}.${SENS}"
LOG="$OUT_DIR/07_make_fixed_markers.${RUN_ID}.log"
exec > >(tee -a "$LOG") 2>&1
echo "=== $(date -u +%FT%TZ) 07_make_fixed_markers.sh $RUN_ID ==="

FIXED_SNPLIST="$OUT_DIR/${RUN_ID}.fixed.snplist"
MANIFEST="$OUT_DIR/${RUN_ID}.marker_manifest.tsv"
GENO_MAF_BIM="$OUT_DIR/${RUN_ID}.geno_maf_filtered.bim"
if [[ -e "$FIXED_SNPLIST" && $OVERWRITE -eq 0 ]]; then
  echo "FATAL: $FIXED_SNPLIST already exists. This is a FROZEN marker set -- pass"
  echo "--overwrite only if you are intentionally rebuilding after fixing a bug in"
  echo "this script itself, not to chase a different SNP count for any ancient sample."
  exit 3
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PARAMS_JSON=$(python3 -c "
import sys, json
sys.path.insert(0, '$SCRIPT_DIR')
from lib_ecotype_v2 import load_config, resolve_marker_params
cfg = load_config('$CONFIG')
print(json.dumps(resolve_marker_params(cfg, '$PANEL', '$SENS')))
")
echo "resolved parameters from config: $PARAMS_JSON"
GENO=$(python3 -c "import json,sys; print(json.loads(sys.argv[1])['geno'])" "$PARAMS_JSON")
MAF=$(python3 -c "import json,sys; print(json.loads(sys.argv[1])['maf'])" "$PARAMS_JSON")
LD_WINDOW_KB=$(python3 -c "import json,sys; print(json.loads(sys.argv[1])['ld_window_kb'])" "$PARAMS_JSON")
LD_R2=$(python3 -c "import json,sys; print(json.loads(sys.argv[1])['ld_r2'])" "$PARAMS_JSON")

RAW_SNPS=$(wc -l < "${BFILE}.bim")
echo "raw_snps (full panel .bim, before reference restriction): $RAW_SNPS"

# --- 1. reference-only, optional capture-bait restriction, TV/ALL split ---
STEP1_EXTRACT="$OUT_DIR/${RUN_ID}.step1_extract.snplist"
if [[ "$LIBTYPE" == "capture" ]]; then
  BAIT_OVERLAP=$(wc -l < "$CAPTURE_SNP_LIST")
  echo "bait_overlap_snps: $BAIT_OVERLAP"
  BASE_LIST="$CAPTURE_SNP_LIST"
else
  BAIT_OVERLAP="NA"
  BASE_LIST=""
fi

python3 -c "
import sys
sys.path.insert(0, '$SCRIPT_DIR')
from lib_ecotype_v2 import is_transversion
base = set()
base_list = '$BASE_LIST'
if base_list:
    with open(base_list) as fh:
        base = set(l.strip() for l in fh if l.strip())
track = '$TRACK'
kept = 0
with open('${BFILE}.bim') as fh, open('$STEP1_EXTRACT', 'w') as out:
    for line in fh:
        parts = line.split()
        chrom, snpid, _cm, pos, a1, a2 = parts[:6]
        if base_list and snpid not in base:
            continue
        if track == 'TV':
            tv = is_transversion(a1, a2)
            if tv is not True:
                continue
        out.write(snpid + chr(10))
        kept += 1
print(f'after_TV_ALL: {kept}', file=sys.stderr)
"
AFTER_TV_ALL=$(wc -l < "$STEP1_EXTRACT")
echo "after_TV_ALL: $AFTER_TV_ALL"
if [[ "$AFTER_TV_ALL" -eq 0 ]]; then
  echo "FATAL: 0 SNPs survive the ${LIBTYPE}/${TRACK} extract step. Stopping -- do not"
  echo "relax TV/ALL or the capture-compatible list to force a nonzero count."
  exit 3
fi

# --- 2. reference samples, site missingness, MAF (no LD yet) ---
GENO_MAF_OUT="$OUT_DIR/${RUN_ID}.geno_maf_filtered"
plink2 --bfile "$BFILE" --keep "$KEEP" --extract "$STEP1_EXTRACT" \
  --geno "$GENO" --maf "$MAF" --make-bed --out "$GENO_MAF_OUT"
AFTER_SITE_MISS_AND_MAF=$(wc -l < "${GENO_MAF_OUT}.bim")
echo "after_site_missingness+MAF (combined, plink2 applies both in one pass): $AFTER_SITE_MISS_AND_MAF"
cp "${GENO_MAF_OUT}.bim" "$GENO_MAF_BIM"
echo "wrote pre-LD intermediate bim (for 08's paperlike_5kb route if panel B): $GENO_MAF_BIM"

REF_N=$(wc -l < "$KEEP")

if [[ "$SENS" == "thinning_only" ]]; then
  echo "SENS=thinning_only: stopping after geno/MAF filter, no LD pruning requested"
  echo "(this run only exists to produce $GENO_MAF_BIM for 08)"
  exit 0
fi

# --- 3. LD prune -- PLINK2 native 2-arg kb syntax, per spec section 5 ---
# spec explicitly forbids the PLINK1.9-style 3-arg "100kb 10 0.2" form.
plink2 --bfile "$GENO_MAF_OUT" --indep-pairwise "${LD_WINDOW_KB}kb" "$LD_R2" \
  --out "$OUT_DIR/${RUN_ID}.ld"
plink2 --bfile "$GENO_MAF_OUT" --extract "$OUT_DIR/${RUN_ID}.ld.prune.in" \
  --make-bed --out "$OUT_DIR/${RUN_ID}.pruned"
AFTER_LD=$(wc -l < "$OUT_DIR/${RUN_ID}.pruned.bim")
echo "after_LD_or_thinning: $AFTER_LD"

cut -f2 "$OUT_DIR/${RUN_ID}.pruned.bim" > "$FIXED_SNPLIST"
MD5=$(md5sum "$FIXED_SNPLIST" | cut -d' ' -f1)

{
  echo -e "panel\tlibrary_type\ttrack\tsensitivity\treference_samples_n\traw_snps\tbait_overlap_snps\tafter_TV_ALL\tafter_site_missingness_and_MAF\tafter_LD_or_thinning\tgeno\tmaf\tld_window_kb\tld_r2\tmd5"
  echo -e "${PANEL}\t${LIBTYPE}\t${TRACK}\t${SENS}\t${REF_N}\t${RAW_SNPS}\t${BAIT_OVERLAP}\t${AFTER_TV_ALL}\t${AFTER_SITE_MISS_AND_MAF}\t${AFTER_LD}\t${GENO}\t${MAF}\t${LD_WINDOW_KB}\t${LD_R2}\t${MD5}"
} > "$MANIFEST"

echo "OK: wrote $FIXED_SNPLIST (n=$AFTER_LD, md5=$MD5) and $MANIFEST"
