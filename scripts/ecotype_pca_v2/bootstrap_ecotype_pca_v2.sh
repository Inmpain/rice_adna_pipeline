#!/usr/bin/env bash
# Download an immutable repository snapshot and prepare the ordered workflow.
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: bootstrap_ecotype_pca_v2.sh --ref FULL_GIT_COMMIT --dest PARENT_DIR

Downloads exactly one immutable Inmpain/rice_adna_pipeline commit.  Refuses
branches/tags and refuses to overwrite an existing installation.
EOF
}

REF=""
DEST=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --ref) REF="$2"; shift 2 ;;
    --dest) DEST="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "FATAL: unknown argument $1" >&2; usage; exit 2 ;;
  esac
done

[[ "$REF" =~ ^[0-9a-f]{40}$ ]] || {
  echo "FATAL: --ref must be a full 40-character lowercase Git commit" >&2
  exit 2
}
[[ -n "$DEST" ]] || { echo "FATAL: --dest is required" >&2; exit 2; }

mkdir -p "$DEST"
DEST=$(cd "$DEST" && pwd)
INSTALL_DIR="$DEST/rice_adna_pipeline-$REF"
[[ ! -e "$INSTALL_DIR" ]] || {
  echo "FATAL: refusing to overwrite existing $INSTALL_DIR" >&2
  exit 3
}

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
ARCHIVE="$TMP/repo.tar.gz"
URL="https://github.com/Inmpain/rice_adna_pipeline/archive/$REF.tar.gz"

echo "downloading immutable commit: $URL"
curl -fL "$URL" -o "$ARCHIVE"
tar -xzf "$ARCHIVE" -C "$TMP"
EXTRACTED="$TMP/rice_adna_pipeline-$REF"
[[ -d "$EXTRACTED/scripts/ecotype_pca_v2/workflow" ]] || {
  echo "FATAL: downloaded archive lacks the v2 workflow directory" >&2
  exit 4
}
mv "$EXTRACTED" "$INSTALL_DIR"
printf '%s\n' "$REF" > "$INSTALL_DIR/.source_revision"

CONTROLLER="$INSTALL_DIR/scripts/ecotype_pca_v2/workflow/ecotype_pca_workflow.py"
python3 "$CONTROLLER" validate-plan

echo "INSTALL_DIR=$INSTALL_DIR"
echo "SOURCE_REVISION=$REF"
echo "NEXT: cd $INSTALL_DIR"
echo "NEXT: python3 scripts/ecotype_pca_v2/workflow/ecotype_pca_workflow.py --state-dir /home/scratch/yinmt202607/gene/results/ecotype_pca_v2/workflow_state status"
