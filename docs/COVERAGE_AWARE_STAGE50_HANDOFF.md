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

```bash
test -s "$CIVAN_UNION_SITES"
test -s "$CIVAN_UNION_SITES_TV"
test -s "$CIVAN_REFERENCE_KEEP"
```

Then submit through the workflow wrapper:

```bash
bash scripts/ecotype_pca_v2/workflow/submit_stage.sh 50 "$STATE"
```

The old `50_civan_fixed_marker_prototype` attempt is retained for audit only;
it is not the coverage-aware result.
