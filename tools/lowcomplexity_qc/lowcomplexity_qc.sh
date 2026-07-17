#!/usr/bin/env bash
set -euo pipefail

# =====================================================================
# lowcomplexity_qc.sh
#
# 通用工具: 检测"基因命中"类分析里, 有多少条reads落在参考基因组的
# 低复杂度/重复序列区域内 —— 这类区域(如GC重复的脯氨酸/丙氨酸编码
# 密码子串联)容易导致短read比对算法(BWA aln等)产生"看似唯一, 实则
# 可比对到基因组多处"的假阳性唯一比对, 需要单独标记, 不能直接当作
# 可信的功能位点证据使用。
#
# 三个子命令, 分别对应质控流程的三个阶段:
#   build-mask       用dustmasker给参考基因组标记低复杂度区域(每个参考
#                     基因组只需要跑一次, 生成的mask可以反复复用)
#   check-hits       把一份"基因命中"TSV跟mask做批量交叉比对, 输出每条
#                     命中是否落在低复杂度区, 并给出汇总统计
#   blast-spotcheck  对某一条具体的read做BLAST抽查, 验证它是否真的
#                     "唯一"比对, 还是能匹配到基因组多处(bwa aln的
#                     "唯一"判定不可全信, 这一步用更穷举的搜索交叉验证)
#
# 用法示例见每个子命令的 --help
# =====================================================================

SCRIPT_NAME="$(basename "${BASH_SOURCE[0]}")"

log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2; }
die() { log "ERROR: $*"; exit 1; }

usage_main() {
  cat <<EOF
用法: $SCRIPT_NAME <command> [选项]

命令:
  build-mask        用dustmasker给参考基因组生成低复杂度区域BED + BLAST索引
  check-hits         检查一份基因命中TSV, 标记哪些命中落在低复杂度区域
  blast-spotcheck    对BAM里某个区域的一条read做BLAST抽查, 验证比对唯一性

查看具体子命令的选项:
  $SCRIPT_NAME build-mask --help
  $SCRIPT_NAME check-hits --help
  $SCRIPT_NAME blast-spotcheck --help
EOF
}

# =====================================================================
# 子命令1: build-mask
# =====================================================================
usage_build_mask() {
  cat <<EOF
用法: $SCRIPT_NAME build-mask --ref REF.fa --out-dir DIR

对参考基因组跑一次dustmasker, 生成:
  DIR/lowcomplexity.sorted.bed   低复杂度区域(标准BED格式, 可反复复用)
  DIR/blastdb/ref.*              BLAST核酸数据库(供blast-spotcheck使用)

每个参考基因组只需要跑一次这一步, 后续所有check-hits/blast-spotcheck
都复用这次的产出, 不需要重新生成。

选项:
  --ref PATH       参考基因组fasta (必填)
  --out-dir PATH   输出目录 (必填)
  -h, --help       显示此帮助
EOF
}

cmd_build_mask() {
  local ref="" out_dir=""

  while (( $# > 0 )); do
    case "$1" in
      --ref) ref="$2"; shift 2 ;;
      --out-dir) out_dir="$2"; shift 2 ;;
      -h|--help) usage_build_mask; exit 0 ;;
      *) die "build-mask: 未知选项 $1" ;;
    esac
  done

  [[ -n "$ref" ]] || die "build-mask: 必须指定 --ref"
  [[ -n "$out_dir" ]] || die "build-mask: 必须指定 --out-dir"
  [[ -f "$ref" ]] || die "参考基因组不存在: $ref"

  command -v dustmasker >/dev/null 2>&1 || die "dustmasker 未找到 (属于NCBI blast+套件)"
  command -v makeblastdb >/dev/null 2>&1 || die "makeblastdb 未找到 (属于NCBI blast+套件)"

  mkdir -p "$out_dir/blastdb"

  log "[1/3] 跑dustmasker标记低复杂度区域..."
  dustmasker -in "$ref" -out "$out_dir/dustmasker_raw.interval" -outfmt interval

  log "[2/3] 转换成标准BED格式..."
  awk '
    /^>/{ chr=substr($1,2); next }
    {
      start = $1 - 1
      if (start < 0) start = 0
      print chr"\t"start"\t"$3
    }
  ' "$out_dir/dustmasker_raw.interval" \
    | sort -k1,1 -k2,2n \
    > "$out_dir/lowcomplexity.sorted.bed"

  n_regions=$(wc -l < "$out_dir/lowcomplexity.sorted.bed")
  log "  低复杂度区域数: $n_regions"

  log "[3/3] 建BLAST数据库(供blast-spotcheck使用)..."
  makeblastdb -in "$ref" -dbtype nucl -out "$out_dir/blastdb/ref" >/dev/null

  log "完成。"
  log "  低复杂度BED: $out_dir/lowcomplexity.sorted.bed"
  log "  BLAST数据库: $out_dir/blastdb/ref"
}

# =====================================================================
# 子命令2: check-hits
# =====================================================================
usage_check_hits() {
  cat <<EOF
用法: $SCRIPT_NAME check-hits --hits HITS.tsv --mask MASK.bed --out OUT.tsv [选项]

把一份"基因命中"TSV跟低复杂度区域mask做批量交叉比对, 用整条read的比对跨度
(而不是单点坐标)判断是否有重叠 —— 只用起始坐标单点检查会系统性漏检
(read起点在低复杂度区外, 但read主体跨入区内的情况)。

输入TSV要求: 带表头, 且包含染色体列、比对起始位置列、read长度列。
不同数据集的列顺序可能不同, 用下面几个 --xxx-col 参数指定具体是第几列
(从1开始计数), 不需要固定列顺序。

选项:
  --hits PATH        输入的命中TSV, 带表头 (必填)
  --mask PATH        build-mask生成的 lowcomplexity.sorted.bed (必填)
  --out PATH         输出的标注结果TSV (必填)
  --chr-col N         染色体列号 (默认: 7)
  --pos-col N         比对起始位置列号, 1-based坐标 (默认: 8)
  --len-col N         read长度列号 (默认: 13)
  --label-cols N,M    用来拼装标签的列号, 逗号分隔, 按顺序用|连接 (默认: 1,6)
  -h, --help          显示此帮助

示例(用gene_hits_with_metadata.tsv的默认列结构):
  $SCRIPT_NAME check-hits \\
      --hits gene_hits_with_metadata.tsv \\
      --mask db/asian_rice_panel_index/lowcomplexity/lowcomplexity.sorted.bed \\
      --out gene_hits_lowcomplexity_check.tsv
EOF
}

cmd_check_hits() {
  local hits="" mask="" out=""
  local chr_col=7 pos_col=8 len_col=13 label_cols="1,6"

  while (( $# > 0 )); do
    case "$1" in
      --hits) hits="$2"; shift 2 ;;
      --mask) mask="$2"; shift 2 ;;
      --out) out="$2"; shift 2 ;;
      --chr-col) chr_col="$2"; shift 2 ;;
      --pos-col) pos_col="$2"; shift 2 ;;
      --len-col) len_col="$2"; shift 2 ;;
      --label-cols) label_cols="$2"; shift 2 ;;
      -h|--help) usage_check_hits; exit 0 ;;
      *) die "check-hits: 未知选项 $1" ;;
    esac
  done

  [[ -n "$hits" ]] || die "check-hits: 必须指定 --hits"
  [[ -n "$mask" ]] || die "check-hits: 必须指定 --mask"
  [[ -n "$out" ]] || die "check-hits: 必须指定 --out"
  [[ -f "$hits" ]] || die "找不到命中文件: $hits"
  [[ -f "$mask" ]] || die "找不到低复杂度mask: $mask (先跑 build-mask)"
  command -v bedtools >/dev/null 2>&1 || die "bedtools 未找到"

  local tmp_spans
  tmp_spans=$(mktemp)

  log "从 $hits 提取命中区间(用整条read跨度, 列: chr=$chr_col pos=$pos_col len=$len_col)..."
  tail -n+2 "$hits" | awk -F'\t' -v chr_c="$chr_col" -v pos_c="$pos_col" \
      -v len_c="$len_col" -v label_c="$label_cols" '
    BEGIN {
      n = split(label_c, cols, ",")
    }
    {
      start = $pos_c - 1
      end   = $pos_c - 1 + $len_c
      label = $(cols[1])
      for (i = 2; i <= n; i++) label = label"|"$(cols[i])
      print $chr_c"\t"start"\t"end"\t"label
    }
  ' > "$tmp_spans"

  local n_input
  n_input=$(wc -l < "$tmp_spans")
  log "  提取到 $n_input 条命中区间"

  log "与低复杂度mask做交叉比对..."
  bedtools intersect -a "$tmp_spans" -b "$mask" -c > "$out"
  rm -f "$tmp_spans"

  local n_total n_low
  n_total=$(wc -l < "$out")
  n_low=$(awk -F'\t' '$5>0' "$out" | wc -l)

  log "完成。输出: $out"
  echo ""
  echo "=== 汇总统计 ==="
  echo "总命中数: $n_total"
  echo "落在低复杂度区域的命中数: $n_low"
  awk -F'\t' -v n="$n_total" -v l="$n_low" 'BEGIN{
    if (n>0) printf "占比: %.1f%% (%d / %d)\n", l/n*100, l, n
    else print "占比: NA (无命中数据)"
  }'
  echo ""
  echo "落在低复杂度区域的具体命中(标签列):"
  awk -F'\t' '$5>0{print "  "$4}' "$out"
}

# =====================================================================
# 子命令3: blast-spotcheck
# =====================================================================
usage_blast_spotcheck() {
  cat <<EOF
用法: $SCRIPT_NAME blast-spotcheck --bam BAM --region chr:start-end --blastdb DB_PREFIX

从BAM里某个区域抽取一条read, 对参考基因组做BLAST搜索, 验证这条read是否
真的"唯一"比对, 还是能以相近分数匹配到基因组上多个位置(BWA aln的"唯一"
判定基于启发式搜索, 不完全可信, 需要更穷举的BLAST交叉验证)。

选项:
  --bam PATH          要抽取read的BAM文件 (必填)
  --region STR        染色体区域, 如 chr07:29616768-29616789 (必填)
  --blastdb PATH      build-mask生成的BLAST数据库前缀(不含.nin/.nsq等后缀), 
                       即 <out_dir>/blastdb/ref (必填)
  --read-index N      该区域内第几条read做检测, 从1开始 (默认: 1)
  --evalue N          BLAST的e-value阈值 (默认: 1, 尽量放宽以捕捉所有候选)
  -h, --help          显示此帮助

输出解读:
  只看到1行结果(1个比对位置)          → 唯一比对, 可信度较高
  看到多行结果, bitscore/pident相近   → 可能是重复序列导致的多重比对歧义,
                                          之前认定"唯一"的BWA比对结果需要
                                          打折扣使用
EOF
}

cmd_blast_spotcheck() {
  local bam="" region="" blastdb="" read_index=1 evalue=1

  while (( $# > 0 )); do
    case "$1" in
      --bam) bam="$2"; shift 2 ;;
      --region) region="$2"; shift 2 ;;
      --blastdb) blastdb="$2"; shift 2 ;;
      --read-index) read_index="$2"; shift 2 ;;
      --evalue) evalue="$2"; shift 2 ;;
      -h|--help) usage_blast_spotcheck; exit 0 ;;
      *) die "blast-spotcheck: 未知选项 $1" ;;
    esac
  done

  [[ -n "$bam" ]] || die "blast-spotcheck: 必须指定 --bam"
  [[ -n "$region" ]] || die "blast-spotcheck: 必须指定 --region"
  [[ -n "$blastdb" ]] || die "blast-spotcheck: 必须指定 --blastdb"
  [[ -f "$bam" ]] || die "找不到BAM: $bam"
  command -v samtools >/dev/null 2>&1 || die "samtools 未找到"
  command -v blastn >/dev/null 2>&1 || die "blastn 未找到"

  local seq read_name
  read_name=$(samtools view "$bam" "$region" | sed -n "${read_index}p" | cut -f1)
  seq=$(samtools view "$bam" "$region" | sed -n "${read_index}p" | cut -f10)

  [[ -n "$seq" ]] || die "在 $region 里找不到第 $read_index 条read (该区域可能没有这么多reads)"

  log "抽取read: $read_name (${#seq}bp), 区域: $region"

  local tmp_fa
  tmp_fa=$(mktemp --suffix=.fa)
  echo -e ">${read_name}\n${seq}" > "$tmp_fa"

  log "跑BLAST(task=blastn-short, evalue<=$evalue)..."
  echo ""
  echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tqstart\tqend\tsstart\tsend\tevalue\tbitscore"
  blastn -query "$tmp_fa" -db "$blastdb" -task blastn-short \
      -outfmt "6 qseqid sseqid pident length mismatch qstart qend sstart send evalue bitscore" \
      -evalue "$evalue"

  rm -f "$tmp_fa"
}

# =====================================================================
# 主入口
# =====================================================================
[[ $# -ge 1 ]] || { usage_main; exit 1; }

command="$1"; shift

case "$command" in
  build-mask)       cmd_build_mask "$@" ;;
  check-hits)       cmd_check_hits "$@" ;;
  blast-spotcheck)  cmd_blast_spotcheck "$@" ;;
  -h|--help)        usage_main; exit 0 ;;
  *) die "未知命令: $command (可用: build-mask, check-hits, blast-spotcheck)" ;;
esac
