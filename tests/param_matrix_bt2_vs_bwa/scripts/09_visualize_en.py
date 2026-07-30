#!/usr/bin/env python3
"""
09_visualize.py
Visualization for the 9-combo extraction x mapping matrix test:
  Figure 1: Heatmap (3 extraction methods x 3 mapping tools, total q30 reads)
  Figure 2: Grouped bar chart (16 samples x 9 combos, q30 reads)
Outputs PDF (vector format, ready for reports/publications)
"""
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns

# ---------- Color definitions ----------
# Heatmap: two-color gradient (light -> dark)
HEATMAP_LOW = "#E0E0FE"   # light lavender, RGB(224,224,254)
HEATMAP_HIGH = "#0000A7"  # dark blue, RGB(0,0,166)
heatmap_cmap = LinearSegmentedColormap.from_list(
    "custom_blue", [HEATMAP_LOW, HEATMAP_HIGH]
)

# Bar plot: 9 of the 10 SCI palette colors (dropped the last neutral
# brownish tone #CAC1B8, kept the 9 with better mutual contrast)
SCI_9_COLORS = [
    "#CAD097",  # 202,208,151
    "#C0CAB7",  # 192,202,183
    "#A7C0A2",  # 167,192,162
    "#F3E7C8",  # 243,231,200
    "#FAE3D6",  # 250,227,214
    "#E4D9E8",  # 228,217,232
    "#D0E9E5",  # 208,233,229
    "#C2D3DF",  # 194,211,223
    "#99B9BE",  # 153,185,190
]

INPUT_TSV = "summary/final_mapping_summary.tsv"
OUT_HEATMAP_PDF = "summary/heatmap_q30_total.pdf"
OUT_BARPLOT_PDF = "summary/barplot_16samples_9combos.pdf"

COMBO_ORDER = [
    'bt2_old_extract__bt2old_map', 'bt2_old_extract__bt2new_map', 'bt2_old_extract__bwa_map',
    'bt2_new_extract__bt2old_map', 'bt2_new_extract__bt2new_map', 'bt2_new_extract__bwa_map',
    'bwa_extract__bt2old_map', 'bwa_extract__bt2new_map', 'bwa_extract__bwa_map',
]
COMBO_LABELS = [c.replace('_extract__', ' -> ').replace('_map', '') for c in COMBO_ORDER]


def main():
    df = pd.read_csv(INPUT_TSV, sep='\t')
    df[['extract_method', 'map_tool']] = df['combo'].str.extract(r'(\w+)_extract__(\w+)_map')

    # ---------- Figure 1: 3x3 heatmap ----------
    pivot = df.groupby(['extract_method', 'map_tool'])['dedup_q30'].sum().unstack()
    pivot = pivot.reindex(index=['bt2_old', 'bt2_new', 'bwa'], columns=['bt2old', 'bwa', 'bt2new'])

    fig1, ax1 = plt.subplots(figsize=(6.5, 5))
    sns.heatmap(
        pivot, annot=True, fmt='.0f', cmap=heatmap_cmap,
        cbar_kws={'label': 'Total q30 reads'}, ax=ax1,
        linewidths=0.5, linecolor='white'
    )
    ax1.set_title('Extraction Method x Mapping Tool\nTotal q30 Reads', fontsize=13)
    ax1.set_ylabel('Extraction Method (Stage 1)', fontsize=11)
    ax1.set_xlabel('Mapping Tool (Stage 2)', fontsize=11)
    plt.tight_layout()
    plt.savefig(OUT_HEATMAP_PDF, format='pdf')
    plt.close(fig1)
    print(f"Figure 1 done: {OUT_HEATMAP_PDF}")

    # ---------- Figure 2: 16 samples x 9 combos bar chart ----------
    pivot2 = df.pivot(index='sample', columns='combo', values='dedup_q30')
    pivot2 = pivot2.reindex(columns=COMBO_ORDER)

    fig2, ax2 = plt.subplots(figsize=(16, 6))
    pivot2.plot(kind='bar', ax=ax2, color=SCI_9_COLORS, width=0.85, edgecolor='none')
    ax2.set_ylabel('q30 reads', fontsize=11)
    ax2.set_xlabel('Sample', fontsize=11)
    ax2.set_title('q30 Reads by Sample and Combo (16 samples x 9 combos)', fontsize=13)
    ax2.legend(COMBO_LABELS, bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8, title='Extraction -> Mapping')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(OUT_BARPLOT_PDF, format='pdf')
    plt.close(fig2)
    print(f"Figure 2 done: {OUT_BARPLOT_PDF}")


if __name__ == "__main__":
    main()
