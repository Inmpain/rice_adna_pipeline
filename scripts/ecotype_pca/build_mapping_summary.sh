#!/usr/bin/env bash
set -euo pipefail

# Combines besthit's per-sample kept_reads (codex/oryza-competitive-mapping
# branch's besthit_summary.tsv, produced by submit_oryza_besthit.sh merge)
# with map_besthit_to_irgsp.sh's per-sample [done] line (mapped/
# duplicates_flagged/mapq30/mapq20) into one table.
#
# This is the "比对情况" slice of docs/ECOTYPE_PCA_EXECUTION_PLAN.md section
# 3's Phase 0 coverage census table -- not the full census yet (no genome-
# wide coverage / panel-intersection columns, those need
# summarize_irgsp_coverage.sh, still unwritten). Just besthit input vs
# mapping output, the slice available right after step ① finishes for all
# 16 samples.
#
# besthit_summary.tsv column order (from oryza_besthit_damage_filter.py's
# module docstring, codex/oryza-competitive-mapping branch):
#   sample  input_reads  reads_with_alignment  reads_with_oryza_hit
#   kept_reads  rejected_nonoryza_better  rejected_no_oryza
#   rejected_low_quality  unclassified_reads
# kept_reads is exactly what feeds map_besthit_to_irgsp.sh as
# <sample>.besthit_oryza.fastq.gz -- confirmed against the smoke-test log
# (LV6000619499: besthit kept_reads and "2445 sequences have been
# processed" in the bwa aln log matched exactly).
#
# Usage: build_mapping_summary.sh <besthit_summary.tsv> <mapirgsp_log_dir> <out.tsv>

BESTHIT_SUMMARY="${1:?usage: build_mapping_summary.sh <besthit_summary.tsv> <mapirgsp_log_dir> <out.tsv>}"
LOG_DIR="${2:?see usage above}"
OUT_TSV="${3:?see usage above}"

[[ -f "$BESTHIT_SUMMARY" ]] || { echo "ERROR: besthit summary not found: $BESTHIT_SUMMARY (run submit_oryza_besthit.sh merge first)" >&2; exit 1; }
[[ -d "$LOG_DIR" ]] || { echo "ERROR: log dir not found: $LOG_DIR" >&2; exit 1; }

TMP_MAP="$(mktemp)"
trap 'rm -f "$TMP_MAP"' EXIT

# Parse every [done] line across all mapirgsp/smoke logs (both share the
# same "[done] <sample>: mapped=... duplicates_flagged=... mapq>=30_nondup=...
# mapq>=20_nondup=..." format) into sample<TAB>mapped<TAB>dup<TAB>q30<TAB>q20
grep -h '^\[done\]' "$LOG_DIR"/*.log \
  | sed -E 's/^\[done\] ([^:]+): mapped=([0-9]+) duplicates_flagged=([0-9]+) mapq>=30_nondup=([0-9]+) mapq>=20_nondup=([0-9]+).*/\1\t\2\t\3\t\4\t\5/' \
  > "$TMP_MAP"

n_map=$(wc -l < "$TMP_MAP")
n_uniq=$(cut -f1 "$TMP_MAP" | sort -u | wc -l)
if [[ "$n_map" -ne "$n_uniq" ]]; then
    echo "ERROR: $LOG_DIR has $n_map [done] lines but only $n_uniq unique samples" >&2
    echo "  -- likely a re-run left old + new logs mixed together (e.g. a" >&2
    echo "  smoke log and a submit log for the same sample). Clean up the" >&2
    echo "  stale log before trusting this table:" >&2
    cut -f1 "$TMP_MAP" | sort | uniq -d >&2
    exit 1
fi

{
    printf 'sample\tbesthit_input_reads\tbesthit_kept_reads\tmapped\tmapped_pct_of_kept\tduplicates_flagged\tdup_pct_of_mapped\tmapq30_nondup\tmapq20_nondup\n'
    tail -n +2 "$BESTHIT_SUMMARY" | while IFS=$'\t' read -r sample input_reads reads_with_alignment reads_with_oryza_hit kept_reads rest; do
        map_row="$(awk -F'\t' -v s="$sample" '$1==s{print; exit}' "$TMP_MAP")"
        if [[ -z "$map_row" ]]; then
            echo "WARNING: no mapping [done] line found for sample $sample -- skipping" >&2
            continue
        fi
        mapped="$(cut -f2 <<< "$map_row")"
        dup="$(cut -f3 <<< "$map_row")"
        q30="$(cut -f4 <<< "$map_row")"
        q20="$(cut -f5 <<< "$map_row")"
        mapped_pct=$(awk -v m="$mapped" -v k="$kept_reads" 'BEGIN{ if (k>0) printf "%.1f", 100*m/k; else print "NA" }')
        dup_pct=$(awk -v d="$dup" -v m="$mapped" 'BEGIN{ if (m>0) printf "%.1f", 100*d/m; else print "NA" }')
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$sample" "$input_reads" "$kept_reads" "$mapped" "$mapped_pct" "$dup" "$dup_pct" "$q30" "$q20"
    done
} > "$OUT_TSV"

echo "[done] wrote $OUT_TSV" >&2
