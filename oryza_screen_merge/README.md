# Merge existing Oryza candidate FASTQs

This Snakemake workflow performs one operation only: for each sample, merge the
existing candidate-Oryza FASTQs from shotgun, panel1, and panel2 into one
compressed FASTQ.

```text
shotgun candidate FASTQ ---+
panel1 candidate FASTQ ----+--> combined candidate-Oryza FASTQ
panel2 candidate FASTQ ----+
```

It does not run BWA, create or read BAM files, remap to IRGSP, remove
duplicates, apply MAPQ30, calculate coverage, call variants, or count gene hits.

The committed configuration matches the current server layout:

```text
/home/scratch/yinmt202607/tests/param_matrix_bt2_vs_bwa/01.reads_combined/bwa/
├── <sample>.prefiltered.IRGSP1.mapped.fq
├── <sample>_RicePanel1.bwa.primary_mapped.fastq.gz
└── <sample>_RicePanel2.bwa.primary_mapped.fastq.gz
```

## Configure

The real input paths and filename patterns are already filled in. Each source
pattern contains a `{sample}` wildcard, so the three files above are grouped by
their common sample ID.

The `{sample}` value extracted from the three patterns must agree. The workflow
stops with a readable error if a source is empty or a sample is missing from any
of the three sources.

The output defaults to:

```text
/home/scratch/yinmt202607/gene/results/oryza_candidates_combined
```

## Download without touching the active repository

The existing `rice_adna_pipeline` checkout contains active untracked work, so do
not clone into it or switch its branch. Download only these workflow files into
the separate gene scripts directory:

```bash
mkdir -p /home/scratch/yinmt202607/gene/scripts/oryza_screen_merge
cd /home/scratch/yinmt202607/gene/scripts/oryza_screen_merge

base_url="https://raw.githubusercontent.com/Inmpain/rice_adna_pipeline/codex/oryza-screen-merge/oryza_screen_merge"

curl -L --retry 8 --retry-connrefused --connect-timeout 30 \
  -o Snakefile "$base_url/Snakefile"
curl -L --retry 8 --retry-connrefused --connect-timeout 30 \
  -o config.yaml "$base_url/config.yaml"
curl -L --retry 8 --retry-connrefused --connect-timeout 30 \
  -o README.md "$base_url/README.md"
```

## Check the planned jobs

```bash
cd /home/scratch/yinmt202607/gene/scripts/oryza_screen_merge
snakemake -n -p
```

## Run

```bash
snakemake --cores 1 --printshellcmds
```

For Slurm, use the cluster's configured Snakemake executor/profile.

## Final output

```text
<output_dir>/combined_mapped_fastq/<sample>.oryza_candidates.combined.fastq.gz
```

The workflow reads the plain shotgun `.fq` directly, decompresses the panel1
and panel2 `.fastq.gz` files, and recompresses all three streams into one gzip
file. This is more portable across downstream tools.
