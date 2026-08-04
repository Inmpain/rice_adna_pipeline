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

## Configure

Edit `config.yaml` and replace the panel1 and panel2 placeholder paths. Each
source pattern must include a `{sample}` wildcard. Source-specific filenames
are supported, for example:

```yaml
panel1:
  directory: /path/to/panel1/candidate_fastq
  pattern: "{sample}_RicePanel1.asian_rice_panel.primary_mapped.fastq.gz"
```

The `{sample}` value extracted from the three patterns must agree. The workflow
stops with a readable error if a source is empty or a sample is missing from any
of the three sources.

## Check the planned jobs

```bash
cd /path/to/rice_adna_pipeline/oryza_screen_merge
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

The workflow decompresses the three inputs and recompresses them as one gzip
stream instead of merely joining compressed byte streams. This is more portable
across downstream tools.
