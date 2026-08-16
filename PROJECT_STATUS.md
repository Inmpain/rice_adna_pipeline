# Rice aDNA pipeline — project status

Last updated: 2026-08-16. This file was previously truncated; this compact
version restores the current project-level view. Workstream details remain in
their branch handoff documents rather than being duplicated here.

## Stable results

- The 3×3 extraction/mapping comparison is complete. BWA at the extraction
  stage is the decisive improvement; see `docs/09_extraction_mapping_matrix_final.md`.
- The ecotype-PCA population-label and UNK-filtering preparation exists for all
  three modern panels. The v1 16×3 plots are explicitly first-look diagnostics,
  not final or cross-sample-comparable PCA results.

## Active workstream: ecotype PCA v2

- Active branch: `codex/ecotype-pca-panel`.
- Statistical design: reference-first, one frozen modern marker set and one PCA
  coordinate system per panel; ancient samples are projection-only.
- Implemented analysis components: scripts 00–08 (input validation, manifests,
  PLINK conversion, audits, reference-set construction, marker freezing and
  Panel-B 5kb thinning).
- Implemented execution control: a receipt-based ordered workflow under
  `scripts/ecotype_pca_v2/workflow/`. It refuses skipped stages, invalidates
  stale receipts after code/config changes, preserves every failed attempt and
  creates returnable debug bundles.
- Corrected on 2026-08-16: exact Civán label strings, EIGENSTRAT REF=2/ALT=0
  specification, TV/ALL seed invariant, and the Panel-B LD forward-halo bug.

## Current next stage

Run stage 00 on the login node, then stage 10 through SLURM; stage 10 scans
large panel metadata and one `samtools flagstat` pass per BAM. The first
server attempt on `33ae004` exposed a PLINK2 synthetic-fixture incompatibility
(10 axis builders, while current PLINK2 requires at least 50 for LD); the
fixture is now 60 and PLINK threads are allocation-bounded. After stage 10,
run stage 20 (full modern Civán sanity PCA) inside SLURM. Stage 30 requires
human review of PC1–2 and PC3–4.

## Fail-closed blockers

- v2 scripts 09–20 are not implemented; fixed-marker ancient prototype and
  production stages remain locked in `workflow.json`.
- Panel B raw 720 versus filtered 718 must be decided from technical evidence,
  not population labels.
- Capture-track work has no confirmed bait BED. Shotgun work is independent.
- Whole-genus versus taxonomic-tier ancient BAM input must be audited before
  production interpretation.

See `docs/ECOTYPE_PCA_PANEL.md` for handoff state,
`docs/ECOTYPE_PCA_V2_WORKFLOW.md` for execution, and
`docs/ECOTYPE_PCA_V2_SPEC.md` for frozen statistical details.
