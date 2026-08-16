# Coverage-aware Civán Stage 50 handoff

## Current scientific state

Stage 00--40 and the original Stage 50 prototype were executed. The original
prototype selected the first 1000 transversion records in file order. Those
markers were concentrated near the beginning of chromosome 1 and are not a
genome-wide marker sample. Its `0/1000` callable result must not be used as a
biological conclusion.

The coverage survey was subsequently run across all 16 ancient BAMs. The
coverage-aware marker universe is based on ancient coverage intersected with
the Civán panel. ALL and TV are separate tracks; TV remains the damage-robust
primary track and ALL is damage-sensitive sensitivity analysis.

All BAMs are interpreted as `pooled_mixed_capture_plus_shotgun`. No capture or
shotgun separation may be inferred. `inputs.capture_bait_bed` remains null and
capture-bait intersection is blocked until a real BED path is supplied.

Every ancient sample is projected when technically possible. Low marker counts
are reported, not silently discarded:

* `formal_validation_candidate`: callable_n >= 200
* `exploratory_projection`: callable_n 50--199
* `descriptive_only`: callable_n < 50

These labels do not change frozen MAPQ/BaseQ, MAF, or pseudo-haplotype rules.

## Workflow implementation

The registered Stage 50 command is:

```text
scripts/ecotype_pca_v2/workflow/runners/50_civan_coverage_aware_projection.sh
```

It requires these exported inputs:

```text
CIVAN_UNION_SITES
CIVAN_UNION_SITES_TV
CIVAN_REFERENCE_KEEP
CIVAN_ANCIENT_SAMPLES
```

It writes independent `ALL/` and `TV/` outputs, keeps modern PCA axes fixed,
and records technical execution separately from scientific projection status.

## Why an old-looking status may reappear

The receipt controller is fail-closed. Changing workflow definitions,
implementation files, or tracked tests changes the implementation digest. The
affected stage becomes `STALE`, and downstream stages are locked until the
stage is rerun in order. This is expected and is not evidence that the earlier
PCA calculations were lost.

Do not run a runner directly on the login node and do not reuse an old attempt
directory. Use the controller/SLURM wrapper so that a new timestamped attempt
and receipt are created.

## Required checks before rerunning Stage 50

Use the one-command server entry. It searches the configured
`results_v2_root` with `find`, accepts only one structurally valid ALL union,
TV union and 595-sample Civán reference keep-list, checks all 16 BAMs, confirms
the `94d5366` QC behavior and coverage-aware workflow registration, reruns the
quick local Stage 40 contract if its receipt is stale, then submits Stage 50.
It never chooses a newest candidate when more than one valid path exists.

```bash
bash scripts/ecotype_pca_v2/workflow/submit_coverage_aware_stage50.sh
```

After the SLURM job finishes it resolves the new attempt from the controller
receipt (never from a hard-coded timestamp), verifies and prints ALL/TV
`pca_qc.tsv`, `scientific_projection.tsv`, and `callability.tsv`, then prints
workflow status. It does not accept or unlock Stage 60.

If discovery reports multiple candidates, inspect the printed real paths,
export only the intended path(s), and rerun the same command. Do not substitute
`/path/to/...` placeholders.

The old `50_civan_fixed_marker_prototype` attempt is retained for audit only;
it is not the coverage-aware result. PCA QC permits an `.ind` sample to be
absent from smartpca's `.evec` only when its label is `Ancient`; omission of a
modern sample remains a hard failure. The omitted ancient IDs are recorded in
`missing_sample_ids`.

Local Claude account/session exports (`users.json`, `login_history.json`, and
unrelated `conversations.json`) are private and must never be uploaded or
committed. Root `.gitignore` rules provide an additional guard.
