#!/usr/bin/env python3
"""
10_best_combo_per_sample.py
针对每个样本单独找出表现最优的组合(避免跨样本对比造成的误读),
分别按 q30 reads 和 gene_hit_reads 两个指标统计, 并汇总每个组合
成为"单样本最优解"的次数, 用于判断BWA+BWA是否是稳定的、普遍的赢家,
而不只是总量堆出来的第一名。
"""
import pandas as pd

INPUT_TSV = "summary/final_mapping_summary.tsv"
OUT_TSV = "summary/best_combo_per_sample.tsv"

df = pd.read_csv(INPUT_TSV, sep="\t")

# ---------- 按 q30 reads 找每个样本的最优组合 ----------
best_q30 = df.loc[df.groupby('sample')['dedup_q30'].idxmax()][
    ['sample', 'combo', 'dedup_q30', 'gene_hit_reads']
].rename(columns={'combo': 'best_combo_by_q30', 'dedup_q30': 'q30_value'})

print("=" * 70)
print("每个样本: q30 reads最优的组合")
print("=" * 70)
print(best_q30.to_string(index=False))

print("\n各组合成为'q30最优'的次数(16个样本里):")
q30_counts = best_q30['best_combo_by_q30'].value_counts()
print(q30_counts.to_string())

# ---------- 按 gene_hit_reads 找每个样本的最优组合 ----------
best_genehit = df.loc[df.groupby('sample')['gene_hit_reads'].idxmax()][
    ['sample', 'combo', 'gene_hit_reads', 'dedup_q30']
].rename(columns={'combo': 'best_combo_by_genehit', 'gene_hit_reads': 'genehit_value'})

print("\n" + "=" * 70)
print("每个样本: gene_hit_reads最优的组合")
print("=" * 70)
print(best_genehit.to_string(index=False))

print("\n各组合成为'gene_hit最优'的次数(16个样本里):")
genehit_counts = best_genehit['best_combo_by_genehit'].value_counts()
print(genehit_counts.to_string())

# ---------- 合并输出成一张表, 存盘 ----------
merged = best_q30.merge(
    best_genehit, on='sample', suffixes=('', '_genehit')
)
merged = merged[['sample', 'best_combo_by_q30', 'q30_value',
                  'best_combo_by_genehit', 'genehit_value']]
merged.to_csv(OUT_TSV, sep="\t", index=False)

print(f"\n完成, 结果已保存: {OUT_TSV}")

# ---------- 额外提示: 两个指标选出的最优组合是否一致 ----------
n_consistent = (merged['best_combo_by_q30'] == merged['best_combo_by_genehit']).sum()
print(f"\nq30最优 与 gene_hit最优 结论一致的样本数: {n_consistent} / {len(merged)}")
