# Oryza candidate-read screen and merge

This Snakemake workflow performs only the candidate-read extraction requested:

1. Screen shotgun, panel1, and panel2 reads independently against
   `asian_rice_panel.fa` using `bwa aln` / `bwa samse`.
2. Keep primary mapped records (`samtools view -F 0x904`).
3. Convert each source BAM back to compressed FASTQ.
4. Concatenate the three candidate FASTQs per sample.

It deliberately does not remap to IRGSP, mark/remove duplicates, apply MAPQ30,
calculate coverage, call variants, or count gene hits.

## Configure

Edit `config.yaml` and replace the reference, panel1, and panel2 placeholder
paths. A source `pattern` may include source-specific text, for example:

```yaml
panel1:
  directory: /path/to/panel1
  pattern: "{sample}_RicePanel1.bbduk.lowcomp_filtered.fq"
```

The workflow intentionally requires the same sample set in all three sources.
It stops with a readable error if any source or sample is missing.

## Check the planned jobs

```bash
cd /path/to/oryza_screen_merge
snakemake -n -p
```

## Run locally

```bash
snakemake --cores 20 --printshellcmds
```

For Slurm, use the cluster's configured Snakemake executor/profile.

## Final output

```text
<output_dir>/combined_mapped_fastq/<sample>.oryza_candidates.combined.fastq.gz
```

The combined FASTQ is rebuilt as one gzip stream rather than merely joining
three compressed byte streams, which is more portable across downstream tools.
