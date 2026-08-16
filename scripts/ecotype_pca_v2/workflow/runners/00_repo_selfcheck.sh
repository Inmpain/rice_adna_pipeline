#!/usr/bin/env bash
# Pure/code-level gate.  No real panel data and no compute-server mutation.
set -euo pipefail

: "${RICE_PCA_REPO_ROOT:?workflow controller must set RICE_PCA_REPO_ROOT}"
: "${RICE_PCA_CONFIG:?workflow controller must set RICE_PCA_CONFIG}"
: "${RICE_PCA_ATTEMPT_DIR:?workflow controller must set RICE_PCA_ATTEMPT_DIR}"

cd "$RICE_PCA_REPO_ROOT"

python3 scripts/ecotype_pca_v2/workflow/ecotype_pca_workflow.py \
  --config "$RICE_PCA_CONFIG" \
  --state-dir "$RICE_PCA_STATE_DIR" \
  validate-plan

for script in scripts/ecotype_pca_v2/*.sh \
              scripts/ecotype_pca_v2/workflow/runners/*.sh \
              scripts/ecotype_pca_v2/bootstrap_ecotype_pca_v2.sh; do
  bash -n "$script"
done

python3 scripts/ecotype_pca_v2/tests/test_lib_ecotype_v2.py
python3 scripts/ecotype_pca_v2/tests/test_06_and_01_logic.py
python3 scripts/ecotype_pca_v2/tests/test_08_streaming.py
python3 scripts/ecotype_pca_v2/tests/test_04_ld_chunking.py
python3 scripts/ecotype_pca_v2/workflow/tests/test_workflow_controller.py
python3 scripts/ecotype_pca_v2/workflow/tests/test_collect_server_evidence.py

python3 - "$RICE_PCA_CONFIG" <<'PY'
import sys
from pathlib import Path

config_path = Path(sys.argv[1])
config = config_path.read_text()
required_config_fragments = [
    "ancient:\n  mapq: 30\n  baseq: 30",
    "panel_A_3k:",
    "  maf: 0.01",
    "panel_C_civan:",
    "    - japonica_(temperate)\n    - japonica_(tropical)",
]
for fragment in required_config_fragments:
    assert fragment in config, f"missing frozen config fragment: {fragment!r}"

spec = Path("docs/ECOTYPE_PCA_V2_SPEC.md").read_text()
assert "REF = 2" in spec and "ALT = 0" in spec
assert "sample + panel + track" not in spec
assert "`track` **不得**进入随机种子" in spec

caller = Path("scripts/ecotype_pca/pseudo_haploid_call.py").read_text()
assert 'if base == ref:' in caller and 'fout.write("2\\n")' in caller
assert 'elif base == alt:' in caller and 'fout.write("0\\n")' in caller
print("PASS: frozen scientific invariants match config/spec/caller")
PY

{
  echo -e "check\tstatus"
  echo -e "plan_schema\tPASS"
  echo -e "shell_syntax\tPASS"
  echo -e "pure_python_tests\tPASS"
  echo -e "frozen_invariants\tPASS"
} > "$RICE_PCA_ATTEMPT_DIR/repo_selfcheck.tsv"

echo "PASS: repository self-check"
