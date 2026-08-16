#!/usr/bin/env bash
set -euo pipefail
ROOT=${RICE_PCA_REPO_ROOT:-$(cd "$(dirname "$0")/../../../.." && pwd)}; cd "$ROOT"
for s in 09_export_fixed_reference_eigenstrat.py 10_call_ancient_fixed_markers.py 11_build_ancient_callability.py 12_build_ancient_overlap_matrix.py 13_merge_ancients_fixed_panel.py 15_pca_qc.py 16_projection_summary.py 17_exact_mask_validation.py 18_validation_metrics.py; do python3 "scripts/ecotype_pca_v2/$s" --help >/dev/null; done
bash scripts/ecotype_pca_v2/14_run_fixed_smartpca.sh --help >/dev/null
python3 -m unittest discover -s scripts/ecotype_pca_v2/tests -p 'test_fixed_projection_contract.py' -v
echo "PASS: 40_fixed_projection_implementation"
