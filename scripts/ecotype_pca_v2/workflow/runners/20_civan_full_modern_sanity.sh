#!/usr/bin/env bash
# Full modern-only Civán PCA.  Must run inside a SLURM allocation.
set -euo pipefail

: "${SLURM_JOB_ID:?this stage must run inside SLURM, not on the login node}"
: "${RICE_PCA_REPO_ROOT:?workflow controller must set RICE_PCA_REPO_ROOT}"
: "${RICE_PCA_CONFIG:?workflow controller must set RICE_PCA_CONFIG}"
: "${RICE_PCA_ATTEMPT_DIR:?workflow controller must set RICE_PCA_ATTEMPT_DIR}"

cd "$RICE_PCA_REPO_ROOT"

eval "$(python3 - "$RICE_PCA_CONFIG" <<'PY'
import shlex
import sys
import yaml

cfg = yaml.safe_load(open(sys.argv[1]))
p = cfg['inputs']['panel_C_civan']
root = p['dir']
prefix = p['prefix']
suffix = p['filtered_suffix']
values = {
    'PANEL_SNP': f'{root}/{prefix}.snp',
    'PANEL_IND': f'{root}/{prefix}{suffix}.ind',
    'PANEL_GENO_ES': f'{root}/{prefix}{suffix}.eigenstratgeno',
    'PANEL_GENO_PLAIN': f'{root}/{prefix}{suffix}.geno',
    'NUM_PCS': cfg['pca']['num_pcs'],
    'NUM_OUTLIER_ITER': cfg['pca']['numoutlieriter'],
    'NUM_CHROM': cfg['pca']['numchrom'],
    'NUM_THREADS': cfg['pca']['numthreads'],
}
for key, value in values.items():
    print(f'{key}={shlex.quote(str(value))}')
PY
)"

[[ -s "$PANEL_SNP" ]] || { echo "FATAL: missing $PANEL_SNP" >&2; exit 2; }
[[ -s "$PANEL_IND" ]] || { echo "FATAL: missing $PANEL_IND" >&2; exit 2; }
if [[ -s "$PANEL_GENO_ES" ]]; then
  PANEL_GENO="$PANEL_GENO_ES"
elif [[ -s "$PANEL_GENO_PLAIN" ]]; then
  PANEL_GENO="$PANEL_GENO_PLAIN"
else
  echo "FATAL: neither filtered genotype file exists: $PANEL_GENO_ES / $PANEL_GENO_PLAIN" >&2
  exit 2
fi

SNP_N=$(wc -l < "$PANEL_SNP")
GENO_N=$(wc -l < "$PANEL_GENO")
IND_N=$(wc -l < "$PANEL_IND")
[[ "$SNP_N" -eq 2365188 ]] || { echo "FATAL: Civán SNP rows=$SNP_N, expected 2365188" >&2; exit 3; }
[[ "$GENO_N" -eq "$SNP_N" ]] || { echo "FATAL: geno rows=$GENO_N != SNP rows=$SNP_N" >&2; exit 3; }
[[ "$IND_N" -eq 1055 ]] || { echo "FATAL: filtered Civán samples=$IND_N, expected 1055" >&2; exit 3; }

POPLIST="$RICE_PCA_ATTEMPT_DIR/civan_domesticated.poplistname.txt"
python3 - "$RICE_PCA_CONFIG" "$PANEL_IND" "$POPLIST" <<'PY'
import collections
import sys
import yaml

cfg = yaml.safe_load(open(sys.argv[1]))
labels = cfg['panel_C_civan']['axis_labels']
expected = cfg['panel_C_civan']['expected_axis_builder_n']
counts = collections.Counter()
with open(sys.argv[2]) as handle:
    for line in handle:
        fields = line.split()
        if len(fields) < 3:
            raise SystemExit(f'FATAL: malformed .ind row: {line!r}')
        counts[fields[2]] += 1
missing = [label for label in labels if counts[label] == 0]
axis_n = sum(counts[label] for label in labels)
if missing or axis_n != expected:
    raise SystemExit(
        f'FATAL: Civán axis-builder gate failed: missing={missing} axis_n={axis_n} '
        f'expected={expected} counts={dict(counts)}'
    )
with open(sys.argv[3], 'w') as out:
    for label in labels:
        out.write(label + '\n')
print(f'PASS: exactly {axis_n} domesticated axis builders; counts=' +
      repr({label: counts[label] for label in labels}))
PY

python3 - "$PANEL_GENO" "$IND_N" <<'PY'
import sys
path, expected = sys.argv[1], int(sys.argv[2])
with open(path) as handle:
    first = handle.readline().rstrip('\r\n')
    if len(first) != expected:
        raise SystemExit(f'FATAL: first genotype row width={len(first)}, expected={expected}')
print(f'PASS: first genotype row width={expected}')
PY

PREFIX="$RICE_PCA_ATTEMPT_DIR/full_civan"
PAR="$PREFIX.par"
cat > "$PAR" <<EOF
genotypename:    $PANEL_GENO
snpname:         $PANEL_SNP
indivname:       $PANEL_IND
evecoutname:     $PREFIX.evec
evaloutname:     $PREFIX.eval
poplistname:     $POPLIST
lsqproject:      YES
numoutevec:      $NUM_PCS
numoutlieriter:  $NUM_OUTLIER_ITER
numchrom:        $NUM_CHROM
numthreads:      $NUM_THREADS
EOF

smartpca -p "$PAR" > "$PREFIX.smartpca.log" 2>&1
[[ -s "$PREFIX.evec" ]] || { echo "FATAL: smartpca produced no $PREFIX.evec" >&2; exit 4; }
[[ -s "$PREFIX.eval" ]] || { echo "FATAL: smartpca produced no $PREFIX.eval" >&2; exit 4; }

EVEC_N=$(awk 'NR>1{n++} END{print n+0}' "$PREFIX.evec")
[[ "$EVEC_N" -eq "$IND_N" ]] || { echo "FATAL: evec samples=$EVEC_N, expected=$IND_N" >&2; exit 4; }

python3 scripts/ecotype_pca/plot_pca_projection.py \
  --evec "full_civan=$PREFIX.evec" --pc-x 1 --pc-y 2 \
  --title "Full Civán modern panel: PC1-PC2; 595 domesticated axis builders" \
  --out "$PREFIX.PC1_PC2.png"
python3 scripts/ecotype_pca/plot_pca_projection.py \
  --evec "full_civan=$PREFIX.evec" --pc-x 3 --pc-y 4 \
  --title "Full Civán modern panel: PC3-PC4; 595 domesticated axis builders" \
  --out "$PREFIX.PC3_PC4.png"

{
  echo -e "metric\tvalue"
  echo -e "panel_snp\t$PANEL_SNP"
  echo -e "panel_geno\t$PANEL_GENO"
  echo -e "panel_ind\t$PANEL_IND"
  echo -e "snp_rows\t$SNP_N"
  echo -e "geno_rows\t$GENO_N"
  echo -e "modern_samples\t$IND_N"
  echo -e "axis_builders\t595"
  echo -e "evec_samples\t$EVEC_N"
  echo -e "slurm_job_id\t$SLURM_JOB_ID"
} > "$PREFIX.summary.tsv"

sha256sum "$PAR" "$PREFIX.evec" "$PREFIX.eval" \
  "$PREFIX.PC1_PC2.png" "$PREFIX.PC3_PC4.png" "$PREFIX.summary.tsv" \
  > "$PREFIX.sha256"

echo "PASS: full Civán modern-only sanity computation"
echo "REVIEW_REQUIRED: $PREFIX.PC1_PC2.png $PREFIX.PC3_PC4.png $PREFIX.smartpca.log $PREFIX.summary.tsv"
