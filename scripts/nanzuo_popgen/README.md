# 南佐水稻 popgen 管线 (nanzuo_popgen)

南佐遗址(1.nanzuo_shotgun / 2.nanzuo_popgen_yancheng / 6.nanzuo_function_yancheng)
三份数据的 Oryza 筛选 → 合并 → IRGSP 定量比对 → 统计 流水线。

## 概览

```
输入三份(统一 key = YWL1-A3487 全名):

  popgen   = 2.nanzuo_popgen_yancheng/*.bbduk.lowcomp_filtered.fq   (未压缩)
  shotgun  = 1.nanzuo_shotgun/*.taxa_cleaned.fq.gz                  (压缩)
  function = 6.nanzuo_function_yancheng/*.bam

  popgen  ──bwa提取(asian_rice_panel)──> fq.gz ┐
  shotgun ──bwa提取(asian_rice_panel)──> fq.gz ┼─合并──> {sample}.combined.fastq.gz
  function ──bam→fq(不跑bwa)──────────> fq.gz ┘
                                                    │
                                                    ▼  比对 irgsp.fa (bwa | bowtie2 -N1 对比)
                             sort → collate → fixmate -m → markdup(只标记不删, 无 -r)
                                                    │
                                                    ▼
                  {sample}.dedup.bam  +  q20/25/30 计数  +  基因组覆盖度
```

## 路径

```
参考(均已建好 bwa 索引):
  asian_rice_panel.fa  /home/scratch/yinmt202607/db/asian_rice_panel_index/asian_rice_panel.fa   (提取用)
  irgsp.fa             /home/scratch/yinmt202607/db/asian_rice_panel_index/irgsp.fa              (最终比对用)
  irgsp_bt2idx         /home/scratch/yinmt202607/db/asian_rice_panel_index/irgsp_bt2idx         (irgsp 的 bowtie2 索引)

输出 /home/scratch/yinmt202607/nanzuo/:
  00.extract/{popgen,shotgun}/{sample}.bwa.primary_mapped.fastq.gz
  00.extract/function/{sample}.bam2fq.fastq.gz
  01.merge/{sample}.combined.fastq.gz
  02.map_irgsp/{bwa,bt2new}/{sample}.dedup.bam(.bai)   # 重复只标记(0x400), 未删除
  02.map_irgsp/{bwa,bt2new}/stats/{sample}.tsv
  03.stats/read_qc.tsv  +  coverage_summary.tsv
```

## 参数

- **提取(阶段①, 比对 asian_rice_panel.fa)**：`bwa aln -l 1024 -n 0.01 -o 2` +
  `samtools view -F 0x904`(只留 primary mapped)。
- **最终比对(阶段②, 比对 irgsp.fa)**，每个样本跑两份对比：
  - `bwa`：`bwa aln -l 1024 -n 0.01 -o 2`
  - `bt2new`：`bowtie2 -k 3 -L 22 -N 1 -i S,1,1.15 --mp 1,1 --rdg 0,1 --rfg 0,1
    --score-min L,0,-0.1 --no-unal`(即 9 格矩阵测试里的“Bowtie2 新参数 -N1”)
- **去重**：`markdup` 不带 `-r`，只打 0x400 标记、不删除。
- 单端数据，用 `bwa samse` / `bowtie2 -U`。

## 运行顺序

```bash
cd /home/scratch/yinmt202607/nanzuo/scripts

# 1. sbatch 提交提取(bwa) + function bam→fq
bash submit_extract.sh 20
squeue -u $USER            # 等全部跑完

# 2. 合并(登录节点, 不 sbatch)
for f in /home/scratch/yinmt202607/nanzuo/00.extract/popgen/*.bwa.primary_mapped.fastq.gz; do
  s=$(basename "$f" .bwa.primary_mapped.fastq.gz)
  bash merge_sample.sh "$s"
done

# 3. sbatch 提交最终比对(bwa + bt2new)
bash submit_map.sh 20
squeue -u $USER            # 等全部跑完

# 4. 汇总统计(登录节点)
bash collect_stats.sh
```

## 统计表列定义

`read_qc.tsv`:

| 列 | 含义 |
|---|---|
| sample | YWL1-AXXXX |
| mapper | bwa / bt2new |
| merged_reads | 合并后 combined.fastq.gz 的 reads 数 |
| primary_mapped | 比对到 irgsp 的 primary mapped reads(去重前) |
| dup_flagged | 被 markdup 标为重复(0x400)的 reads 数 |
| dup_rate_pct | dup_flagged / primary_mapped × 100 |
| q20 / q25 / q30 | 排除重复(0x400) 且 MAPQ≥20/25/30 的 reads 数(只计数, 不删) |

`coverage_summary.tsv`(基于去重后 reads, 排除 0x400):

| 列 | 含义 |
|---|---|
| cov_bases | 基因组被覆盖的碱基总数 |
| cov_pct | 覆盖碱基 / 基因组总长 × 100 |
| mean_depth | 平均测序深度 |

## 注意事项

- 工具(bwa/bowtie2/samtools)需在 sbatch 任务内可用；脚本已做 `module load`
  兜底，若走 conda/mamba 环境，请改 submit 脚本里 `sbatch` 前加 `--export` 或
  在单脚本头部把 env 的 bin 加进 PATH。
- 队列分区默认 `comp`，资源默认(提取 32G/9h、比对 24G/12h、bam2fq
  8G/2h)，用 `PARTITION=... bash submit_extract.sh` 覆盖分区，其余在 submit
  脚本里直接改。
- 提取阶段 BAM 只作中间产物用完即删，最终只留 FASTQ；重复不删、q 过滤不删，
  全部只出统计表。
