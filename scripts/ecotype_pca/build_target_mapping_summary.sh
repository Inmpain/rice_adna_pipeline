#!/usr/bin/env bash
# Join taxonomic-tier counts to mapping logs for the narrowed target read set.
set -euo pipefail

TIER_SUMMARY="${1:?usage: build_target_mapping_summary.sh <taxonomic_tiers_summary.tsv> <map_log_dir> <out.tsv>}"
LOG_DIR="${2:?see usage above}"
OUT_TSV="${3:?see usage above}"

[[ -f "$TIER_SUMMARY" ]] || { echo "ERROR: missing tier summary: $TIER_SUMMARY" >&2; exit 1; }
[[ -d "$LOG_DIR" ]] || { echo "ERROR: missing log directory: $LOG_DIR" >&2; exit 1; }

mkdir -p "$(dirname "$OUT_TSV")"
tmp_map="$(mktemp "${OUT_TSV}.map.tmp.XXXXXX")"
tmp_out="${OUT_TSV}.tmp.$$"
trap 'rm -f "$tmp_map" "$tmp_out"' EXIT

shopt -s nullglob
logs=("$LOG_DIR"/*.log)
shopt -u nullglob
[[ "${#logs[@]}" -gt 0 ]] || { echo "ERROR: no *.log files under $LOG_DIR" >&2; exit 1; }

grep -h '^\[done\]' "${logs[@]}" \
  | sed -E 's/^\[done\] ([^:]+): mapped=([0-9]+) duplicates_flagged=([0-9]+) mapq>=30_nondup=([0-9]+) mapq>=20_nondup=([0-9]+).*/\1\t\2\t\3\t\4\t\5/' \
  > "$tmp_map"

n_map=$(wc -l < "$tmp_map")
n_uniq=$(cut -f1 "$tmp_map" | sort -u | wc -l)
if [[ "$n_map" -ne "$n_uniq" ]]; then
    echo "ERROR: duplicate [done] lines found for these samples:" >&2
    cut -f1 "$tmp_map" | sort | uniq -d >&2
    exit 1
fi

awk -F'\t' -v OFS='\t' '
    NR==FNR {
        mapped[$1]=$2; dup[$1]=$3; q30[$1]=$4; q20[$1]=$5
        next
    }
    FNR==1 {
        print "sample", "besthit_kept_reads", "target_reads", "other_oryza_reads", \
              "target_pct_of_kept", "mapped", "mapped_pct_of_target", \
              "duplicates_flagged", "dup_pct_of_mapped", "mapq30_nondup", "mapq20_nondup"
        next
    }
    {
        sample=$1
        if (!(sample in mapped)) {
            print "WARNING: no mapping [done] line for " sample > "/dev/stderr"
            next
        }
        mapped_pct=($5 > 0 ? sprintf("%.1f", 100*mapped[sample]/$5) : "NA")
        dup_pct=(mapped[sample] > 0 ? sprintf("%.1f", 100*dup[sample]/mapped[sample]) : "NA")
        print sample, $4, $5, $6, $7, mapped[sample], mapped_pct, \
              dup[sample], dup_pct, q30[sample], q20[sample]
    }
' "$tmp_map" "$TIER_SUMMARY" > "$tmp_out"

mv "$tmp_out" "$OUT_TSV"
trap - EXIT
rm -f "$tmp_map"
echo "[done] wrote $OUT_TSV" >&2
