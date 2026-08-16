#!/usr/bin/env bash
# Real-server, read-only input/tool audit plus synthetic integration tests.
set -euo pipefail

: "${SLURM_JOB_ID:?stage 10 scans large panel files and BAM flags; run it through SLURM, not on a login node}"
: "${RICE_PCA_REPO_ROOT:?workflow controller must set RICE_PCA_REPO_ROOT}"
: "${RICE_PCA_CONFIG:?workflow controller must set RICE_PCA_CONFIG}"
: "${RICE_PCA_ATTEMPT_DIR:?workflow controller must set RICE_PCA_ATTEMPT_DIR}"

export RICE_PCA_PLINK_THREADS="${SLURM_CPUS_PER_TASK:-1}"

cd "$RICE_PCA_REPO_ROOT"
V2="scripts/ecotype_pca_v2"

python3 "$V2/00_validate_inputs.py" \
  --config "$RICE_PCA_CONFIG" \
  --out-dir "$RICE_PCA_ATTEMPT_DIR/input_validation" \
  --allow-capture-blocked

python3 "$V2/01_make_panel_manifest.py" \
  --config "$RICE_PCA_CONFIG" \
  --out-dir "$RICE_PCA_ATTEMPT_DIR/panel_manifest"

bash "$V2/tests/test_integration_synthetic.sh" "$RICE_PCA_REPO_ROOT/$V2" \
  2>&1 | tee "$RICE_PCA_ATTEMPT_DIR/integration_synthetic.log"

python3 scripts/ecotype_pca/tests/test_pseudo_haploid_call.py \
  2>&1 | tee "$RICE_PCA_ATTEMPT_DIR/pseudo_haploid_regression.log"

python3 "$V2/workflow/collect_server_evidence.py" \
  --config "$RICE_PCA_CONFIG" \
  --out-dir "$RICE_PCA_ATTEMPT_DIR/server_evidence"

echo "PASS: server preflight and synthetic integration"
