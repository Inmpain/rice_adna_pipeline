#!/usr/bin/env bash
# Discover the real coverage-aware Civán inputs, validate them, submit a new
# controller-owned Stage 50 attempt, and print the three scientific outputs.
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: submit_coverage_aware_stage50.sh [STATE_DIR]

Default STATE_DIR:
  /home/scratch/yinmt202607/gene/results/ecotype_pca_v2/workflow_state

The script honors already-exported CIVAN_UNION_SITES,
CIVAN_UNION_SITES_TV, and CIVAN_REFERENCE_KEEP.  Otherwise it searches the
versioned config's results_v2_root and proceeds only when exactly one valid
candidate exists for each input.  It never selects a newest file or reuses an
attempt directory when discovery is ambiguous.
EOF
}

case "${1:-}" in
  -h|--help) usage; exit 0 ;;
esac
[[ $# -le 1 ]] || { usage >&2; exit 2; }

die() {
  echo "FATAL: $*" >&2
  exit 2
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
CONFIG="$REPO_ROOT/scripts/ecotype_pca_v2/config/ecotype_pca_v2.yaml"
CTL="$SCRIPT_DIR/ecotype_pca_workflow.py"
STATE_DIR="${1:-/home/scratch/yinmt202607/gene/results/ecotype_pca_v2/workflow_state}"

[[ -f "$CONFIG" ]] || die "config not found: $CONFIG"
[[ -f "$CTL" ]] || die "controller not found: $CTL"
[[ -f "$SCRIPT_DIR/workflow.json" ]] || die "workflow plan not found"
mkdir -p "$STATE_DIR"
STATE_DIR="$(cd "$STATE_DIR" && pwd)"
cd "$REPO_ROOT"

CONFIG_EXPORTS="$(python3 - "$CONFIG" <<'PY'
import shlex
import sys

import yaml

cfg = yaml.safe_load(open(sys.argv[1]))
if cfg['inputs'].get('capture_bait_bed') is not None:
    raise SystemExit('FATAL: capture_bait_bed must remain null for this pooled analysis')
panel = cfg['inputs']['panel_C_civan']
values = {
    'RESULTS_V2_ROOT': cfg['results_v2_root'],
    'ANCIENT_BAM_DIR': cfg['inputs']['ancient_bam_dir'],
    'PANEL_IND': f"{panel['dir']}/{panel['prefix']}{panel['filtered_suffix']}.ind",
    'EXPECTED_REFERENCE_N': cfg['panel_C_civan']['expected_axis_builder_n'],
}
for key, value in values.items():
    print(f'{key}={shlex.quote(str(value))}')
PY
)" || exit $?
eval "$CONFIG_EXPORTS"

[[ -d "$RESULTS_V2_ROOT" ]] || die "results_v2_root does not exist: $RESULTS_V2_ROOT"
[[ -s "$PANEL_IND" ]] || die "Civán filtered .ind is missing or empty: $PANEL_IND"

valid_sites_file() {
  python3 - "$1" <<'PY' >/dev/null 2>&1
import csv
import sys

required = {'snp_id', 'chrom', 'pos', 'n_samples_covered', 'samples_covered'}
with open(sys.argv[1], newline='') as handle:
    reader = csv.DictReader(handle, delimiter='\t')
    if not reader.fieldnames or not required.issubset(reader.fieldnames):
        raise SystemExit(1)
    ids = []
    for row in reader:
        ids.append(row['snp_id'])
        int(row['pos'])
        int(row['n_samples_covered'])
if not ids or len(ids) != len(set(ids)):
    raise SystemExit(1)
PY
}

valid_reference_keep() {
  python3 - "$1" "$PANEL_IND" "$EXPECTED_REFERENCE_N" <<'PY' >/dev/null 2>&1
import sys

keep_path, ind_path, expected_text = sys.argv[1:]
expected = int(expected_text)
keep = []
with open(keep_path) as handle:
    for line in handle:
        fields = line.split()
        if fields:
            keep.append(fields[1] if len(fields) >= 2 else fields[0])
ind_ids = {line.split()[0] for line in open(ind_path) if line.split()}
if len(keep) != expected or len(keep) != len(set(keep)):
    raise SystemExit(1)
if not set(keep).issubset(ind_ids):
    raise SystemExit(1)
PY
}

discover_unique() {
  local variable="$1"
  local label="$2"
  local pattern="$3"
  local validator="$4"
  local supplied="${!variable:-}"
  local candidate
  local candidates=()

  if [[ -n "$supplied" ]]; then
    [[ "$supplied" != *'/path/to/'* ]] || die "$variable still contains a placeholder"
    [[ "$supplied" == /* ]] || die "$variable must be an absolute real path: $supplied"
    [[ -s "$supplied" ]] || die "$variable is missing or empty: $supplied"
    "$validator" "$supplied" || die "$variable failed structural validation: $supplied"
    printf -v "$variable" '%s' "$supplied"
    export "$variable"
    return
  fi

  while IFS= read -r candidate; do
    [[ -s "$candidate" ]] || continue
    if "$validator" "$candidate"; then
      candidates+=("$candidate")
    fi
  done < <(find "$RESULTS_V2_ROOT" -type f -name "$pattern" -print 2>/dev/null | LC_ALL=C sort)

  if [[ ${#candidates[@]} -eq 0 ]]; then
    die "no valid $label found below $RESULTS_V2_ROOT"
  fi
  if [[ ${#candidates[@]} -ne 1 ]]; then
    echo "FATAL: found ${#candidates[@]} valid candidates for $label; refusing to guess:" >&2
    printf '  %s\n' "${candidates[@]}" >&2
    echo "Export $variable to the intended real path, then rerun this same script." >&2
    exit 2
  fi

  printf -v "$variable" '%s' "${candidates[0]}"
  export "$variable"
}

discover_unique CIVAN_UNION_SITES "ALL union sites" "ancient_union_sites.tsv" valid_sites_file
discover_unique CIVAN_UNION_SITES_TV "TV union sites" "ancient_union_sites.TV.tsv" valid_sites_file
discover_unique CIVAN_REFERENCE_KEEP "595-sample Civán reference keep-list" "*.reference_samples.keep" valid_reference_keep

export CIVAN_ANCIENT_SAMPLES="${CIVAN_ANCIENT_SAMPLES:-LV6000619499 LV6000619917 LV6000620016 LV6000620032 LV6000620166 LV6000620172 LV6000654686 LV6000654698 LV7008416272 LV7008416280 LV7008416294 LV7008416329 LV7008416339 LV7008416349 LV7008416379 LV7008416407}"

python3 - "$CIVAN_UNION_SITES" "$CIVAN_UNION_SITES_TV" "$CIVAN_ANCIENT_SAMPLES" <<'PY'
import csv
import sys

all_path, tv_path, sample_text = sys.argv[1:]
expected_samples = set(sample_text.split())
if len(expected_samples) != 16:
    raise SystemExit(f'FATAL: expected 16 unique ancient sample IDs, found {len(expected_samples)}')

def load(path):
    with open(path, newline='') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    ids = {row['snp_id'] for row in rows}
    covered = set()
    for row in rows:
        covered.update(value for value in row['samples_covered'].split(',') if value)
    return ids, covered

all_ids, all_samples = load(all_path)
tv_ids, tv_samples = load(tv_path)
if not tv_ids.issubset(all_ids):
    raise SystemExit('FATAL: TV union marker IDs are not a subset of ALL union marker IDs')
if all_samples != expected_samples:
    raise SystemExit('FATAL: ALL union samples do not match the requested 16 ancient samples: '
                     f'missing={sorted(expected_samples-all_samples)} extra={sorted(all_samples-expected_samples)}')
if tv_samples != expected_samples:
    raise SystemExit('FATAL: TV union samples do not match the requested 16 ancient samples: '
                     f'missing={sorted(expected_samples-tv_samples)} extra={sorted(tv_samples-expected_samples)}')
print(f'PASS: coverage inputs agree (ALL={len(all_ids)} markers; TV={len(tv_ids)} markers; samples=16)')
PY

sample_n=0
seen_samples=" "
for sample in $CIVAN_ANCIENT_SAMPLES; do
  [[ "$sample" =~ ^LV[0-9]+$ ]] || die "invalid ancient sample ID: $sample"
  [[ "$seen_samples" != *" $sample "* ]] || die "duplicate ancient sample ID: $sample"
  seen_samples+="$sample "
  sample_n=$((sample_n + 1))
  [[ -s "$ANCIENT_BAM_DIR/$sample.besthit_oryza.irgsp.bam" ]] || \
    die "ancient BAM is missing or empty: $ANCIENT_BAM_DIR/$sample.besthit_oryza.irgsp.bam"
done
[[ $sample_n -eq 16 ]] || die "expected 16 ancient samples, found $sample_n"

grep -q 'missing_sample_ids' scripts/ecotype_pca_v2/15_pca_qc.py || \
  die "15_pca_qc.py predates commit 94d5366's missing-sample fix"
python3 - "$SCRIPT_DIR/workflow.json" <<'PY'
import json
import sys

plan = json.load(open(sys.argv[1]))
stage = next(item for item in plan['stages'] if item['id'] == '50_civan_fixed_marker_prototype')
expected = ['bash', 'scripts/ecotype_pca_v2/workflow/runners/50_civan_coverage_aware_projection.sh']
if stage['command'] != expected:
    raise SystemExit(f"FATAL: Stage 50 is not registered to the coverage-aware runner: {stage['command']}")
print('PASS: Stage 50 is registered to the coverage-aware runner')
PY

echo "STATE_DIR=$STATE_DIR"
echo "CIVAN_UNION_SITES=$CIVAN_UNION_SITES"
echo "CIVAN_UNION_SITES_TV=$CIVAN_UNION_SITES_TV"
echo "CIVAN_REFERENCE_KEEP=$CIVAN_REFERENCE_KEEP"
echo "CIVAN_ANCIENT_SAMPLES=$CIVAN_ANCIENT_SAMPLES"
test -s "$CIVAN_UNION_SITES"
test -s "$CIVAN_UNION_SITES_TV"
test -s "$CIVAN_REFERENCE_KEEP"

python3 "$CTL" --state-dir "$STATE_DIR" status
NEXT_OUTPUT="$(python3 "$CTL" --state-dir "$STATE_DIR" next)" || {
  printf '%s\n' "$NEXT_OUTPUT" >&2
  exit 2
}
printf '%s\n' "$NEXT_OUTPUT"
NEXT_ID="$(printf '%s\n' "$NEXT_OUTPUT" | sed -n 's/^NEXT=\([^ ]*\).*/\1/p')"

if [[ "$NEXT_ID" == "40_fixed_projection_implementation" ]]; then
  echo "Stage 40 is stale because its tracked QC implementation changed; rerunning its local tests."
  python3 "$CTL" --state-dir "$STATE_DIR" run 40_fixed_projection_implementation
  NEXT_OUTPUT="$(python3 "$CTL" --state-dir "$STATE_DIR" next)" || {
    printf '%s\n' "$NEXT_OUTPUT" >&2
    exit 2
  }
  printf '%s\n' "$NEXT_OUTPUT"
  NEXT_ID="$(printf '%s\n' "$NEXT_OUTPUT" | sed -n 's/^NEXT=\([^ ]*\).*/\1/p')"
fi

[[ "$NEXT_ID" == "50_civan_fixed_marker_prototype" ]] || \
  die "controller says the next stage is ${NEXT_ID:-unknown}, not Stage 50; no job was submitted"

# Slurm's documented default is --export=ALL; set the equivalent environment
# option explicitly so all four validated CIVAN_* values reach the job.
export SBATCH_EXPORT=ALL
bash "$SCRIPT_DIR/submit_stage.sh" 50 "$STATE_DIR"

RECEIPT="$STATE_DIR/receipts/50_civan_fixed_marker_prototype.json"
[[ -s "$RECEIPT" ]] || die "Stage 50 returned without a receipt: $RECEIPT"
ATTEMPT_DIR="$(python3 - "$RECEIPT" <<'PY'
import json
import sys
print(json.load(open(sys.argv[1]))['attempt_dir'])
PY
)"
[[ -d "$ATTEMPT_DIR" ]] || die "receipt attempt directory not found: $ATTEMPT_DIR"

python3 - "$ATTEMPT_DIR" "$CIVAN_ANCIENT_SAMPLES" <<'PY'
import csv
import sys
from pathlib import Path

attempt = Path(sys.argv[1])
expected_samples = set(sys.argv[2].split())
for track in ('ALL', 'TV'):
    root = attempt / track
    paths = {
        'pca_qc': root / f'civan.{track}.pca_qc.tsv',
        'scientific_projection': root / f'civan.{track}.scientific_projection.tsv',
        'callability': root / f'civan.{track}.callability.tsv',
    }
    for label, path in paths.items():
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f'FATAL: missing {track} {label}: {path}')

    with paths['pca_qc'].open(newline='') as handle:
        qc = {row['metric']: row['value'] for row in csv.DictReader(handle, delimiter='\t')}
    missing = {value for value in qc['missing_sample_ids'].split(',') if value}
    if not missing.issubset(expected_samples):
        raise SystemExit(f'FATAL: {track} pca_qc lists non-ancient missing IDs: {sorted(missing-expected_samples)}')
    if int(qc['evec_samples_missing_from_pca']) != len(missing):
        raise SystemExit(f'FATAL: {track} pca_qc missing count/ID list disagree')

    with paths['callability'].open(newline='') as handle:
        callability = list(csv.DictReader(handle, delimiter='\t'))
    with paths['scientific_projection'].open(newline='') as handle:
        scientific = list(csv.DictReader(handle, delimiter='\t'))
    if {row['sample'] for row in callability} != expected_samples:
        raise SystemExit(f'FATAL: {track} callability does not contain exactly the 16 requested samples')
    if {row['sample'] for row in scientific} != expected_samples:
        raise SystemExit(f'FATAL: {track} scientific_projection does not contain exactly the 16 requested samples')
    for row in scientific:
        n = int(row['callable_n'])
        expected = ('formal_validation_candidate' if n >= 200 else
                    'exploratory_projection' if n >= 50 else 'descriptive_only')
        if row['technical_execution'] != 'PASS' or row['scientific_projection'] != expected:
            raise SystemExit(f"FATAL: {track}/{row['sample']} classification mismatch")
print(f'PASS: verified ALL/TV pca_qc, scientific_projection, and callability for all 16 samples')
PY

echo "NEW_ATTEMPT=$ATTEMPT_DIR"
for track in ALL TV; do
  echo "===== $track pca_qc.tsv ====="
  sed -n '1,20p' "$ATTEMPT_DIR/$track/civan.$track.pca_qc.tsv"
  echo "===== $track scientific_projection.tsv ====="
  sed -n '1,25p' "$ATTEMPT_DIR/$track/civan.$track.scientific_projection.tsv"
  echo "===== $track callability.tsv ====="
  sed -n '1,25p' "$ATTEMPT_DIR/$track/civan.$track.callability.tsv"
done

python3 "$CTL" --state-dir "$STATE_DIR" status
echo "PASS: coverage-aware Stage 50 complete. Stage 60 was not accepted or unlocked by this script."
