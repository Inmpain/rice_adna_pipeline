#!/usr/bin/env bash
set -euo pipefail
if [[ ${1:-} == --help || ${1:-} == -h ]]; then
  echo "usage: $0 --evec-label LABEL --evec PATH --title-prefix TEXT --out-prefix PREFIX"
  echo "  writes PREFIX.PC1_PC2.png .. PREFIX.PC9_PC10.png (5 pairs) via scripts/ecotype_pca/plot_pca_projection.py"
  exit 0
fi

LABEL=""; EVEC=""; TITLE_PREFIX=""; OUT_PREFIX=""
REPO_ROOT="${RICE_PCA_REPO_ROOT:-.}"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --evec-label) LABEL="$2"; shift 2 ;;
    --evec) EVEC="$2"; shift 2 ;;
    --title-prefix) TITLE_PREFIX="$2"; shift 2 ;;
    --out-prefix) OUT_PREFIX="$2"; shift 2 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done
for v in LABEL EVEC TITLE_PREFIX OUT_PREFIX; do
  [[ -n "${!v}" ]] || { echo "missing required argument (see --help)" >&2; exit 2; }
done
[[ -s "$EVEC" ]] || { echo "FATAL: evec not found or empty: $EVEC" >&2; exit 3; }

for PAIR in "1 2" "3 4" "5 6" "7 8" "9 10"; do
  read -r X Y <<<"$PAIR"
  python3 "$REPO_ROOT/scripts/ecotype_pca/plot_pca_projection.py" \
    --evec "${LABEL}=${EVEC}" \
    --pc-x "$X" --pc-y "$Y" \
    --title "${TITLE_PREFIX}: PC${X}-PC${Y}" \
    --out "${OUT_PREFIX}.PC${X}_PC${Y}.png"
done
echo "PASS: wrote 5 PC-pair plots with prefix ${OUT_PREFIX}"
