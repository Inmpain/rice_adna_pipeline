#!/usr/bin/env bash
# 重建 all_wgs_asian_irgsp.acc2taxid（2026-08-08第三次修正：
# 1) 确认 asian_rice_panel.acc2taxid 里 4529/4530 是整份文件系统性反标，
#    不是只有np7等12个基因组——凡是标4529(rufipogon)的应为4530(sativa)，
#    凡是标4530的应为4529，全文件做一次二元互换即可，不需要按基因组名字
#    分别判断。
# 2) 服务器dry-run实测发现 irgsp.acc2taxid 有一行表头
#    (accession/accession.version/taxid的字面表头行)，12行真实数据全部
#    正确标成4530(sativa)，taxid方向本身没问题，但表头行如果被原样cat
#    进最终的1100万行合并文件，第3列会是字面字符串"taxid"而不是数字，
#    besthit脚本按taxid解析这一列时会出问题——merge阶段现在会过滤掉任何
#    第3列不是纯数字的行（对三个源文件都生效，不只是irgsp）。
#
# 三个源文件：
#   1. WGS真核库taxid（独立大文件，406M，跟panel目录不在一起）：
#      /home/database/ref20250728/taxonomy_CPH/wgs_eukaryota.acc2taxid
#   2. 亚洲水稻panel taxid：<panel目录>/asian_rice_panel.acc2taxid
#      （本脚本会先把这份文件里的4529/4530做整体互换，见下方Step 1）
#   3. IRGSP taxid：<panel目录>/irgsp.acc2taxid
#      （本脚本只做统计报告，不自动改——是否也有同样的反标问题还没验证过，
#      见下方Step 2的输出，需要人工确认后再决定要不要修）
#
# 背景：docs/asian_rice_panel_reference_design_conversation.md +
#       docs/ORYZA_BESTHIT_HANDOFF.md 第0.5节
#
# 用法：
#   bash rebuild_all_wgs_asian_irgsp_acc2taxid.sh \
#     --wgs /home/database/ref20250728/taxonomy_CPH/wgs_eukaryota.acc2taxid \
#     --panel-dir /home/scratch/yinmt202607/db/asian_rice_panel_index \
#     [--apply]
# 默认dry-run，只统计不写文件；确认无误后加 --apply 才真正落盘替换。

set -euo pipefail

WGS=""
PANEL_DIR=""
APPLY=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --wgs) WGS="$2"; shift 2 ;;
    --panel-dir) PANEL_DIR="$2"; shift 2 ;;
    --apply) APPLY=1; shift ;;
    *) echo "unknown arg: $1" >&2; exit 1 ;;
  esac
done

if [[ -z "$WGS" || -z "$PANEL_DIR" ]]; then
  echo "usage: $0 --wgs <wgs_eukaryota.acc2taxid路径> --panel-dir <asian_rice_panel_index目录> [--apply]" >&2
  exit 1
fi

PANEL_ACC2TAXID="$PANEL_DIR/asian_rice_panel.acc2taxid"
IRGSP_ACC2TAXID="$PANEL_DIR/irgsp.acc2taxid"
MERGED_ACC2TAXID="$PANEL_DIR/all_wgs_asian_irgsp.acc2taxid"

for f in "$WGS" "$PANEL_ACC2TAXID" "$IRGSP_ACC2TAXID"; do
  [[ -f "$f" ]] || { echo "ERROR: 找不到 $f" >&2; exit 1; }
done

TS="$(date +%Y%m%d-%H%M%S)"

echo "== Step 1: 整体互换 asian_rice_panel.acc2taxid 里的 4529/4530 =="
echo "-- 修正前 taxid(第3列) 分布 --"
cut -f3 "$PANEL_ACC2TAXID" | sort | uniq -c

FIXED_PANEL="$PANEL_DIR/asian_rice_panel.acc2taxid.fixed"
awk 'BEGIN{FS=OFS="\t"}
{
    if ($3 == 4529) { $3 = 4530; changed++ }
    else if ($3 == 4530) { $3 = 4529; changed++ }
    print
}
END { print "changed_lines=" (changed+0) > "/dev/stderr" }' "$PANEL_ACC2TAXID" > "$FIXED_PANEL"

echo "-- 修正后 taxid(第3列) 分布 --"
cut -f3 "$FIXED_PANEL" | sort | uniq -c

echo "-- 抽样对比(前5行，格式: accession  修正前taxid -> 修正后taxid) --"
paste <(cut -f1 "$PANEL_ACC2TAXID" | head -5) \
      <(cut -f3 "$PANEL_ACC2TAXID" | head -5) \
      <(cut -f3 "$FIXED_PANEL" | head -5) | \
  awk 'BEGIN{FS=OFS="\t"} {print $1, $2, "->", $3}'

echo ""
echo "== Step 2: 检查 irgsp.acc2taxid 的taxid分布(不自动改，只报告) =="
echo "-- irgsp.acc2taxid 前5行 --"
head -5 "$IRGSP_ACC2TAXID"
echo "-- irgsp.acc2taxid 第3列(taxid)分布 --"
cut -f3 "$IRGSP_ACC2TAXID" | sort | uniq -c
echo "⚠️ 结合asian_rice_panel.acc2taxid的教训——如果上面显示的物种和taxid对不上"
echo "   （比如整批应该是4530(sativa)的却标成4529，或者反过来），说明"
echo "   irgsp.acc2taxid 也有同样的反标问题，需要另外处理——这个脚本目前"
echo "   不会自动修正它，请把上面的输出贴回来确认后再决定"

echo ""
echo "== Step 2b: 检查三个源文件里第3列不是纯数字的行(比如表头行) =="
echo "这些行如果原样并入合并文件，taxid列会变成非数字字符串，besthit脚本"
echo "解析taxid时会出问题——merge阶段(--apply)会自动过滤掉这些行，不会写入"
echo "最终的 all_wgs_asian_irgsp.acc2taxid"
for label_file in "wgs_eukaryota.acc2taxid:$WGS" "asian_rice_panel.acc2taxid(修正后):$FIXED_PANEL" "irgsp.acc2taxid:$IRGSP_ACC2TAXID"; do
  label="${label_file%%:*}"
  f="${label_file#*:}"
  bad_count=$(awk -F'\t' '$3 !~ /^[0-9]+$/' "$f" | wc -l | tr -d ' ')
  echo "-- ${label} 非数字taxid行数 = ${bad_count} --"
  if [[ "$bad_count" -gt 0 ]]; then
    awk -F'\t' '$3 !~ /^[0-9]+$/' "$f" | head -5
  fi
done

echo ""
echo "== Step 3: 三源合并统计(dry-run先看行数，不代表最终结果) =="
echo "wgs_eukaryota.acc2taxid 行数: $(wc -l < "$WGS")"
echo "asian_rice_panel.acc2taxid(修正后) 行数: $(wc -l < "$FIXED_PANEL")"
echo "irgsp.acc2taxid 行数: $(wc -l < "$IRGSP_ACC2TAXID")"
if [[ -f "$MERGED_ACC2TAXID" ]]; then
  echo "现有 all_wgs_asian_irgsp.acc2taxid 行数(即将被替换): $(wc -l < "$MERGED_ACC2TAXID")"
fi

if [[ "$APPLY" -eq 0 ]]; then
  echo ""
  echo "== DRY RUN 完成，没有写任何文件 =="
  echo "确认上面的数字、taxid互换结果、irgsp.acc2taxid的taxid分布都没问题后，"
  echo "加 --apply 重跑一次才会真正落盘替换"
  rm -f "$FIXED_PANEL"
  exit 0
fi

echo ""
echo "== Step 4 (--apply): 备份原文件，重建合并文件 =="
[[ -f "$MERGED_ACC2TAXID" ]] && cp "$MERGED_ACC2TAXID" "$MERGED_ACC2TAXID.bak-$TS"
cp "$PANEL_ACC2TAXID" "$PANEL_ACC2TAXID.bak-$TS"

mv "$FIXED_PANEL" "$PANEL_ACC2TAXID"

REBUILT="$MERGED_ACC2TAXID.rebuilt"
# 只保留第3列是纯数字的行(丢弃表头等格式不对的行，见Step 2b)
awk -F'\t' '$3 ~ /^[0-9]+$/ {print}' "$WGS" "$PANEL_ACC2TAXID" "$IRGSP_ACC2TAXID" > "$REBUILT"
mv "$REBUILT" "$MERGED_ACC2TAXID"

echo "已完成。备份文件："
[[ -f "$MERGED_ACC2TAXID.bak-$TS" ]] && echo "  $MERGED_ACC2TAXID.bak-$TS"
echo "  $PANEL_ACC2TAXID.bak-$TS"
echo "新 all_wgs_asian_irgsp.acc2taxid 行数: $(wc -l < "$MERGED_ACC2TAXID")"
