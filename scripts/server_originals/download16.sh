#!/usr/bin/env bash
set -euo pipefail

OUT_DIR="${1:-asian_rice_panel_download}"
INCLUDE="${INCLUDE:-genome,gff3}"

ACCESSIONS=(
  GCF_001433935.1
  GCA_009831315.1
  GCA_009830595.1
  GCA_009831275.1
  GCA_009831255.1
  GCA_001623345.2
  GCA_009914875.1
  GCA_009831045.1
  GCA_009831025.1
  GCA_009831355.1
  GCA_009829395.1
  GCA_009831295.1
  GCA_009829375.1
  GCA_001623365.2
  GCA_001952365.3
  GCA_009831335.1
)

usage() {
  cat <<EOF
Usage:
  bash $(basename "$0") [out_dir]

Environment variables:
  INCLUDE   Files to download, default: genome,gff3

Example:
  bash $(basename "$0")
  INCLUDE=genome bash $(basename "$0") rice_panel_genome_only
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

if ! command -v datasets >/dev/null 2>&1; then
  echo "Error: datasets command not found in PATH." >&2
  exit 1
fi

mkdir -p "$OUT_DIR"/zips "$OUT_DIR"/unzipped "$OUT_DIR"/logs
FAILED_FILE="$OUT_DIR/failed_accessions.txt"
SUCCESS_FILE="$OUT_DIR/success_accessions.txt"
: > "$FAILED_FILE"
: > "$SUCCESS_FILE"

total="${#ACCESSIONS[@]}"
i=0

for acc in "${ACCESSIONS[@]}"; do
  i=$((i + 1))
  zip_file="$OUT_DIR/zips/${acc}.zip"
  extract_dir="$OUT_DIR/unzipped/${acc}"
  log_file="$OUT_DIR/logs/${acc}.log"

  echo "[$i/$total] $acc"

  if [[ -s "$zip_file" && -d "$extract_dir/ncbi_dataset" ]]; then
    echo "  skip existing"
    printf '%s\n' "$acc" >> "$SUCCESS_FILE"
    continue
  fi

  if datasets download genome accession \
    "$acc" \
    --include "$INCLUDE" \
    --filename "$zip_file" \
    >"$log_file" 2>&1; then
    mkdir -p "$extract_dir"
    unzip -qo "$zip_file" -d "$extract_dir" >>"$log_file" 2>&1
    printf '%s\n' "$acc" >> "$SUCCESS_FILE"
    echo "  done"
  else
    printf '%s\n' "$acc" >> "$FAILED_FILE"
    echo "  failed, see $log_file" >&2
  fi
done

echo
echo "Finished."
echo "Success list: $SUCCESS_FILE"
echo "Failed list:  $FAILED_FILE"

