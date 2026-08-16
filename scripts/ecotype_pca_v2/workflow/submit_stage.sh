#!/usr/bin/env bash
# Submit guarded compute stages with versioned SLURM resources.
set -euo pipefail

usage() {
  echo "Usage: $0 {10|20|50} [STATE_DIR]"
  echo "Default STATE_DIR: /home/scratch/yinmt202607/gene/results/ecotype_pca_v2/workflow_state"
}

case "${1:-}" in
  -h|--help) usage; exit 0 ;;
esac

STAGE="${1:-}"
STATE_DIR="${2:-/home/scratch/yinmt202607/gene/results/ecotype_pca_v2/workflow_state}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
CTL="$SCRIPT_DIR/ecotype_pca_workflow.py"

command -v sbatch >/dev/null || { echo "FATAL: sbatch is not on PATH" >&2; exit 1; }
mkdir -p "$STATE_DIR"
STATE_DIR="$(cd "$STATE_DIR" && pwd)"

case "$STAGE" in
  10|10_server_preflight)
    STAGE_ID="10_server_preflight"
    JOB_NAME="pca2_preflight"
    MEMORY="8G"
    ;;
  20|20_civan_full_modern_sanity)
    STAGE_ID="20_civan_full_modern_sanity"
    JOB_NAME="pca2_civan_full"
    MEMORY="24G"
    ;;
  50|50_civan_fixed_marker_prototype)
    STAGE_ID="50_civan_fixed_marker_prototype"
    JOB_NAME="pca2_civan_proto"
    MEMORY="8G"
    ;;
  *)
    echo "FATAL: stage must be 10, 20, or 50" >&2
    usage
    exit 2
    ;;
esac

WRAP="$(printf 'cd %q && python3 %q --state-dir %q run %q' \
  "$REPO_ROOT" "$CTL" "$STATE_DIR" "$STAGE_ID")"

echo "submitting $STAGE_ID from $REPO_ROOT"
set +e
sbatch --wait -p comp --exclude=node05,node06 \
  -c 2 --mem "$MEMORY" -t 24:00:00 \
  -J "$JOB_NAME" \
  -o "$STATE_DIR/${STAGE_ID}.%j.slurm.log" \
  --wrap "$WRAP"
SBATCH_RC=$?
set -e

LATEST_LOG=""
for log in "$STATE_DIR/${STAGE_ID}."*.slurm.log; do
  [[ -f "$log" ]] || continue
  if [[ -z "$LATEST_LOG" || "$log" -nt "$LATEST_LOG" ]]; then
    LATEST_LOG="$log"
  fi
done
if [[ -n "$LATEST_LOG" ]]; then
  echo "LATEST_LOG=$LATEST_LOG"
  tail -n 120 "$LATEST_LOG"
fi

python3 "$CTL" --state-dir "$STATE_DIR" status
exit "$SBATCH_RC"
