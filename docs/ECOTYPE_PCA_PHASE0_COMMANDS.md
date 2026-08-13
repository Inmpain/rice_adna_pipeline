# Ecotype PCA Phase 0: ORSC split, remapping, and panel-overlap census

This runbook starts from the completed genus-wide best-hit outputs. It does not
rerun competitive mapping and does not overwrite the existing whole-Oryza BAMs.

The default target is the panel-relevant ORSC set: *O. rufipogon* (4529),
*O. sativa* (4530), and *O. nivara* (4536). This is intentionally narrower than
"all AA-genome Oryza".

## 1. Download the branch scripts

```bash
cd /home/scratch/yinmt202607/gene/scripts
BASE=https://raw.githubusercontent.com/Inmpain/rice_adna_pipeline/codex/ecotype-pca-panel

curl -fsSLO "$BASE/scripts/oryza_besthit/split_besthit_taxonomic_tiers.py"
curl -fsSLO "$BASE/scripts/oryza_besthit/run_split_taxonomic_tiers.sh"
curl -fsSLO "$BASE/scripts/ecotype_pca/map_besthit_to_irgsp.sh"
curl -fsSLO "$BASE/scripts/ecotype_pca/summarize_panel_overlap.py"
curl -fsSLO "$BASE/scripts/ecotype_pca/build_target_mapping_summary.sh"
chmod +x split_besthit_taxonomic_tiers.py run_split_taxonomic_tiers.sh \
  map_besthit_to_irgsp.sh summarize_panel_overlap.py build_target_mapping_summary.sh
```

## 2. Split all completed best-hit samples

```bash
./run_split_taxonomic_tiers.sh check
./run_split_taxonomic_tiers.sh run all
./run_split_taxonomic_tiers.sh merge

column -t -s $'\t' \
  /home/scratch/yinmt202607/gene/results/oryza_competitive_mapping/taxonomic_tiers/taxonomic_tiers_summary.tsv
```

Stop and inspect the by-species table before remapping, especially for
LV7008416349 and LV7008416379:

```bash
column -t -s $'\t' \
  /home/scratch/yinmt202607/gene/results/oryza_competitive_mapping/taxonomic_tiers/taxonomic_tiers_by_species.tsv \
  | less -S
```

## 3. Map the target reads to IRGSP in a new directory

```bash
export BESTHIT_DIR=/home/scratch/yinmt202607/gene/results/oryza_competitive_mapping/taxonomic_tiers
export INPUT_SUFFIX=.target_orsc.fastq.gz
export READSET_LABEL=target_orsc
export OUT_DIR=/home/scratch/yinmt202607/gene/results/ecotype_pca/bam_irgsp_orsc

./map_besthit_to_irgsp.sh check
./map_besthit_to_irgsp.sh smoke LV7008416379
# Inspect the smoke log/output, then:
./map_besthit_to_irgsp.sh submit all
```

After jobs finish:

```bash
./build_target_mapping_summary.sh \
  /home/scratch/yinmt202607/gene/results/oryza_competitive_mapping/taxonomic_tiers/taxonomic_tiers_summary.tsv \
  /home/scratch/yinmt202607/gene/results/ecotype_pca/bam_irgsp_orsc/logs \
  /home/scratch/yinmt202607/gene/results/ecotype_pca/bam_irgsp_orsc/mapping_summary.tsv
```

## 4. Measure real panel overlap on two smoke samples

The output reports Q0 and Q20 side by side, always with BQ20 and duplicate /
secondary / supplementary / QC-fail reads excluded.

```bash
module load python/ 2>/dev/null || true

for sample in LV7008416379 LV7008416294; do
  python3 summarize_panel_overlap.py \
    --sample "$sample" \
    --bam "/home/scratch/yinmt202607/gene/results/ecotype_pca/bam_irgsp_orsc/${sample}.target_orsc.irgsp.bam" \
    --panel civan=/home/scratch/yinmt202607/db/paper1/civan_snp.snp \
    --panel 3k_29m=/home/scratch/yinmt202607/db/29M_3k/NB_final_snp.snp \
    --panel wild720=/home/scratch/yinmt202607/db/6.7M_720/asn720.6m.snp \
    --min-baseq 20 --low-mapq 0 --high-mapq 20 \
    --out "/home/scratch/yinmt202607/gene/results/ecotype_pca/bam_irgsp_orsc/${sample}.panel_overlap.tsv"
done
```

`callable_all_low/high` count covered positions with at least one base matching
either panel allele. `callable_tv_low/high` additionally exclude transitions.
These overlap counts determine whether a PCA smoke test is informative; read
counts alone do not.
