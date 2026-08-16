# `workflow/` — ordered ecotype PCA v2 execution

Purpose: turn the branch handoff order into machine-enforced stages with
content-addressed receipts and reproducible debug attempts.

Contents:

- `workflow.json`: authoritative linear stage/gate list and tracked files.
- `ecotype_pca_workflow.py`: standard-library controller (`status`, `next`,
  `run`, `accept`, `debug-bundle`).
- `runners/`: exact commands for stages currently implemented.
- `collect_server_evidence.py`: read-only raw/filtered sample diff and BAM pair
  flag census used to resolve the 720/718 and overlap questions.
- `tests/`: controller ordering/digest regression tests.

Attempt products and receipts do not belong in the repository. They live in
the user-supplied server `--state-dir`, normally
`/home/scratch/yinmt202607/gene/results/ecotype_pca_v2/workflow_state`.

See `docs/ECOTYPE_PCA_V2_WORKFLOW.md` for complete operator commands and
`docs/ECOTYPE_PCA_V2_SPEC.md` for statistical details.
