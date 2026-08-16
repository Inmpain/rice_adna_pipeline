# Rice aDNA pipeline — project status

Last updated: 2026-08-17. This file was previously truncated; this compact
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
- Implemented analysis components: scripts 00–22, including fixed-reference
  export, ancient calling/callability, merge/projection/QC, coverage survey,
  transversion filtering and scientific projection tiers.
- Implemented execution control: a receipt-based ordered workflow under
  `scripts/ecotype_pca_v2/workflow/`. It refuses skipped stages, invalidates
  stale receipts after code/config changes, preserves every failed attempt and
  creates returnable debug bundles.
- Corrected on 2026-08-16: exact Civán label strings, EIGENSTRAT REF=2/ALT=0
  specification, TV/ALL seed invariant, and the Panel-B LD forward-halo bug.

## Current next stage

Stages 00–40 are complete. Rerun coverage-aware Stage 50 through
`workflow/submit_coverage_aware_stage50.sh`; it discovers and validates the
real ALL union, TV union, Civán reference keep-list and 16 BAMs, then creates a
new controller attempt. The old first-1000-TV-marker attempt remains an audit
artifact and cannot support biological conclusions. TV is the primary track;
ALL is damage-sensitive sensitivity analysis. Do not unlock Stage 60.

## Fail-closed blockers

- Stage 60 remains blocked pending accepted Civán results and exact-mask
  encoding/orientation validation; a technically successful Stage 50 does not
  auto-accept it.
- Panel B raw 720 versus filtered 718 must be decided from technical evidence,
  not population labels.
- Capture-track work has no confirmed bait BED. Shotgun work is independent.
- Whole-genus versus taxonomic-tier ancient BAM input must be audited before
  production interpretation.

See `docs/ECOTYPE_PCA_PANEL.md` for handoff state,
`docs/ECOTYPE_PCA_V2_WORKFLOW.md` for execution, and
`docs/ECOTYPE_PCA_V2_SPEC.md` for frozen statistical details.
