#!/usr/bin/env bash
# 修正 asian_rice_panel.acc2taxid 里的 taxid 错标问题，并把修正结果补丁式
# 合并回 all_wgs_asian_irgsp.acc2taxid（不重建整份文件，只替换这12个基因组
# 对应的行，因为原始三源合并的脚本没有留存，不知道WGS/IRGSP各自那部分是
# 怎么拼的——只动我们确认有问题的这一块，风险最小）。
#
# 背景（详见 docs/asian_rice_panel_reference_design_conversation.md）：
# 用染色体长度+逐条染色体MD5比对确认 asian_rice_panel.fa 里的 np7.Chr1-12
# 与 irgsp.fa 的 chr01-12 完全一致(IDENTICAL)，即 np7 = Nipponbare/IRGSP-1.0
# = O. sativa，但 asian_rice_panel.acc2taxid 里 np7.* 全部被标成了4529
# (O. rufipogon，普通野生稻)——taxid反了。
#
# ⚠️只有 np7 做过这种逐染色体MD5的独立验证。其余11个基因组(mh63/X24_kas/
# azu/arc/liuxu/7个G25_ruf_W*)的物种归属是根据文件命名和来源判断的、
# 沿用对话里的结论，没有再逐一做MD5验证——如果之后发现某个命名也对不上，
# 这个脚本需要跟着改。
#
# 用法：
#   bash fix_asian_rice_panel_taxid.sh /path/to/asian_rice_panel_index目录
# 默认按 dry-run 模式跑，只打印会改多少行、改成什么，不写文件；
# 确认无误后加 --apply 才真正落盘。

set -euo pipefail

DIR="${1:-}"
APPLY=0
for arg in "$@"; do
  [[ "$arg" == "--apply" ]] && APPLY=1
done

if [[ -z "$DIR" || ! -d "$DIR" ]]; then
  echo "usage: $0 <asian_rice_panel_index目录路径> [--apply]" >&2
  exit 1
fi

PANEL_ACC2TAXID="$DIR/asian_rice_panel.acc2taxid"
MERGED_ACC2TAXID="$DIR/all_wgs_asian_irgsp.acc2taxid"

for f in "$PANEL_ACC2TAXID" "$MERGED_ACC2TAXID"; do
  [[ -f "$f" ]] || { echo "ERROR: 找不到 $f" >&2; exit 1; }
done

TS="$(date +%Y%m%d-%H%M%S)"

echo "== Step 1: 生成修正后的 asian_rice_panel.acc2taxid =="
FIXED_PANEL="$DIR/asian_rice_panel.acc2taxid.fixed"
awk 'BEGIN{FS=OFS="\t"}
{
    if ($1 ~ /^(np7|mh63|X24_kas|azu|arc|liuxu)\./) {
        if ($3 != 4530) { changed++ }
        $3 = 4530
    }
    else if ($1 ~ /^G25_ruf_(W1214|W0169|W1750|W3037|W1536|W1726|W2064)\./) {
        if ($3 != 4529) { changed++ }
        $3 = 4529
    }
    print
}
END { print "changed_lines=" (changed+0) > "/dev/stderr" }' "$PANEL_ACC2TAXID" > "$FIXED_PANEL"

echo "-- 修正前后对比(只显示np7的，其余11个基因组同理) --"
grep -E '^np7\.' "$PANEL_ACC2TAXID" | head -3
echo "  ↓ 修正为 ↓"
grep -E '^np7\.' "$FIXED_PANEL" | head -3

echo ""
echo "== Step 2: 把这12个基因组对应的行，在 all_wgs_asian_irgsp.acc2taxid 里做替换 =="
PREFIX_REGEX='^(np7|mh63|X24_kas|azu|arc|liuxu|G25_ruf_(W1214|W0169|W1750|W3037|W1536|W1726|W2064))\.'
OLD_COUNT=$(grep -cE "$PREFIX_REGEX" "$MERGED_ACC2TAXID" || true)
NEW_COUNT=$(wc -l < "$FIXED_PANEL")
echo "旧merged文件里这12个基因组的行数: $OLD_COUNT"
echo "修正后panel文件的总行数(应替换成这个数): $NEW_COUNT"

PATCHED_MERGED="$DIR/all_wgs_asian_irgsp.acc2taxid.patched"
grep -vE "$PREFIX_REGEX" "$MERGED_ACC2TAXID" > "$PATCHED_MERGED"
cat "$FIXED_PANEL" >> "$PATCHED_MERGED"

TOTAL_OLD=$(wc -l < "$MERGED_ACC2TAXID")
TOTAL_NEW=$(wc -l < "$PATCHED_MERGED")
echo "合并文件总行数: 修正前=$TOTAL_OLD 修正后=$TOTAL_NEW (应该相等，如果不等说明有问题，不要往下apply)"

if [[ "$APPLY" -eq 0 ]]; then
  echo ""
  echo "== DRY RUN 完成，没有写任何文件 =="
  echo "结果预览文件："
  echo "  $FIXED_PANEL"
  echo "  $PATCHED_MERGED"
  echo "确认上面的数字和np7示例都对得上以后，加 --apply 重跑一次才会真正落盘"
  exit 0
fi

echo ""
echo "== Step 3 (--apply): 备份原文件并落盘 =="
cp "$PANEL_ACC2TAXID" "$PANEL_ACC2TAXID.bak-$TS"
cp "$MERGED_ACC2TAXID" "$MERGED_ACC2TAXID.bak-$TS"
mv "$FIXED_PANEL" "$PANEL_ACC2TAXID"
mv "$PATCHED_MERGED" "$MERGED_ACC2TAXID"
echo "已完成。备份文件："
echo "  $PANEL_ACC2TAXID.bak-$TS"
echo "  $MERGED_ACC2TAXID.bak-$TS"
