# Ecotype PCA v2 ordered workflow

This is the functional overview for executing and debugging PCA v2. Exact
statistical parameters and acceptance definitions remain in
`ECOTYPE_PCA_V2_SPEC.md`; they are deliberately not repeated as editable shell
arguments here.

## Why this controller exists

The earlier runbooks are valuable history but can be executed out of order, and
chat memory can accidentally mix v1 results, proposed commands and real server
results. `scripts/ecotype_pca_v2/workflow/ecotype_pca_workflow.py` replaces
that execution surface with a version-controlled, fail-closed state machine.

It guarantees the following mechanical properties:

1. Only the first incomplete stage can run. A stage cannot be skipped by
   naming a later ID.
2. Success requires exit code 0 or an explicit evidence-backed manual PASS.
3. A receipt is tied to the stage definition, every tracked implementation
   file, the frozen config and upstream receipt hashes.
4. Changing config/code makes the affected receipt stale and locks downstream
   work. Re-downloading a corrected commit does not silently bless an old run.
5. Each retry uses a new timestamped attempt directory. Failed output is never
   overwritten and is automatically packaged as `*.debug.tar.gz`.
6. Scientific visual review is a manual gate. The controller never decides
   that clusters “look right.”

This controls order and provenance; it cannot prove that an external program or
input dataset is scientifically correct. That is why preflight, exact-mask
validation and manual gates remain separate stages.

## Fixed stage order

| Stage | Purpose | Initial availability |
|---|---|---|
| 00 | Repository/config/spec/pure-test self-check | open |
| 10 | Real server preflight, synthetic integration, BAM/panel evidence in SLURM | open after 00 |
| 20 | Full 2.365M-marker Civán modern-only PCA in SLURM | open after 00/10 |
| 30 | Human review of PC1–2, PC3–4, log and summary | manual gate |
| 40 | Implement/test v2 scripts 09–18 | complete on current server state; rerun if digest is stale |
| 50 | Coverage-aware Civán ALL/TV projection of all 16 ancient samples | open after 40 |
| 60 | 3K MAF=0.01 fixed-marker prototype | locked |
| 70 | Decide 720/718 and run corrected Panel-B audit | locked |
| 80 | All-ancient production plus sensitivity/cross-panel summary | locked |

The machine-readable authority is `workflow/workflow.json`. Later commits open
stages 40 onward only when their implementation and tests actually exist.

## First server run

Use an immutable 40-character commit supplied with the GitHub handoff. The
bootstrap refuses branches and refuses to overwrite an installation:

```bash
bash bootstrap_ecotype_pca_v2.sh \
  --ref <FULL_COMMIT> \
  --dest /home/scratch/yinmt202607/gene/workflow_sources
```

Then define reusable paths:

```bash
SRC=/home/scratch/yinmt202607/gene/workflow_sources/rice_adna_pipeline-<FULL_COMMIT>
STATE=/home/scratch/yinmt202607/gene/results/ecotype_pca_v2/workflow_state
CTL="$SRC/scripts/ecotype_pca_v2/workflow/ecotype_pca_workflow.py"
cd "$SRC"
```

Inspect, then run stage 00 on the login node:

```bash
python3 "$CTL" --state-dir "$STATE" status
python3 "$CTL" --state-dir "$STATE" next
python3 "$CTL" --state-dir "$STATE" run 00_repo_selfcheck
```

Do not run the next stage directly on the login node. Stage 10 reads the 29M
panel metadata and scans BAM flags, so both the controller and runner require a
SLURM allocation. It permits capture to remain explicitly BLOCKED but requires
the shotgun inputs, labels, tools, synthetic integration and pseudo-haploid
regression test to pass.

## Stages 10 and 20 through SLURM

Run the preflight with the versioned submission helper. It requests 2 CPUs,
8 GB memory and 24 hours, waits for completion, then prints workflow status:

```bash
bash scripts/ecotype_pca_v2/workflow/submit_stage.sh 10
```

The helper prints the newest SLURM log tail even when the job fails, so the
`DEBUG_BUNDLE=...` line is visible without another long shell command.

The first real attempt on commit `33ae004` correctly failed before BAM scanning:
PLINK2 2.0 rejected LD estimation on the old 10-axis-builder synthetic fixture.
The corrected fixture has 60 axis builders; production commands do not use
PLINK2's unsafe `--bad-ld` override. BAM paired/proper-pair evidence now uses
one `samtools flagstat` scan per BAM instead of three full scans.

After stage 10 has a valid receipt, stage 20 uses the same helper with 2 CPUs,
24 GB memory and 24 hours:

```bash
bash scripts/ecotype_pca_v2/workflow/submit_stage.sh 20
```

The computational receipt is written only after smartpca, row-count checks and
both PNGs finish successfully. It does not pass the scientific review gate.

## Manual review gate

After reviewing the two PNGs and log, record an evidence-backed decision. Do
not accept merely because the job exited zero:

```bash
ATTEMPT=$(find "$STATE/attempts/20_civan_full_modern_sanity" -mindepth 1 -maxdepth 1 -type d | sort | tail -1)
python3 "$CTL" --state-dir "$STATE" accept 30_civan_full_modern_review \
  --decision PASS \
  --evidence "$ATTEMPT/full_civan.PC1_PC2.png" \
  --evidence "$ATTEMPT/full_civan.PC3_PC4.png" \
  --evidence "$ATTEMPT/full_civan.smartpca.log" \
  --evidence "$ATTEMPT/full_civan.summary.tsv" \
  --note "Reviewed against Civán modern population structure; <write concrete observations>."
```

If the result is not credible, do not create a PASS receipt. Return the debug
bundle and images for code/data review.

## Debug loop

Any failed command prints a line like:

```text
DEBUG_BUNDLE=/.../20260816T....10_server_preflight.debug.tar.gz
```

Send that tarball back. It contains only workflow metadata, config, receipts,
logs and small attempt outputs—not BAMs or genotype matrices. After a fix is
committed, install the new immutable commit but reuse the same `STATE` path.
The digest mechanism will show exactly which receipt became `STALE`; unchanged
accepted stages remain complete.

To recreate a bundle manually for the current stage:

```bash
python3 "$CTL" --state-dir "$STATE" debug-bundle
```

Never edit a receipt by hand and never copy a receipt between state roots.

## Coverage-aware Stage 50

The historical first-1000-TV-marker attempt is diagnostic only. The current
registered Stage 50 fixes modern Civán axes, projects every requested ancient
sample on separate ALL and TV callable marker sets, and records scientific
tiers without changing frozen MAPQ/BaseQ, MAF, or low-information thresholds.
TV is primary; ALL is a damage-sensitive sensitivity track.

From an immutable source snapshot at the current commit, run:

```bash
bash scripts/ecotype_pca_v2/workflow/submit_coverage_aware_stage50.sh
```

The entry script discovers real inputs under the configured results root and
refuses zero or ambiguous candidates. It validates the three input files, all
16 BAMs, the current QC fix and controller registration before submission.
Every retry is created by the controller in a new timestamped attempt. The
script checks and prints both tracks' `pca_qc.tsv`,
`scientific_projection.tsv`, and `callability.tsv` after success. It never
accepts or unlocks Stage 60.
