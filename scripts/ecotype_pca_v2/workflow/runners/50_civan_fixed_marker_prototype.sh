#!/usr/bin/env bash
set -euo pipefail
: "${SLURM_JOB_ID:?stage 50 must run in SLURM}"
: "${RICE_PCA_REPO_ROOT:?}"; : "${RICE_PCA_CONFIG:?}"; : "${RICE_PCA_ATTEMPT_DIR:?}"
cd "$RICE_PCA_REPO_ROOT"
eval "$(python3 - "$RICE_PCA_CONFIG" <<'PY'
import shlex,sys,yaml
c=yaml.safe_load(open(sys.argv[1])); p=c['inputs']['panel_C_civan']
for k,v in {'SNP':f"{p['dir']}/{p['prefix']}{p['filtered_suffix']}.snp",'BAMDIR':c['inputs']['ancient_bam_dir']}.items(): print(f'{k}={shlex.quote(v)}')
PY
)"
BAM="$BAMDIR/LV7008416379.besthit_oryza.irgsp.bam"; [[ -s "$SNP" && -s "$BAM" ]] || { echo "FATAL: missing Civán SNP or prototype BAM" >&2; exit 2; }
FIXED="$RICE_PCA_ATTEMPT_DIR/civan.TV.prototype.snp"
python3 - "$SNP" "$FIXED" <<'PY'
import sys
tv={frozenset(('A','C')),frozenset(('A','T')),frozenset(('C','G')),frozenset(('G','T'))}; n=0
with open(sys.argv[1]) as src,open(sys.argv[2],'w') as out:
 for line in src:
  f=line.split()
  if len(f)==6 and frozenset((f[4].upper(),f[5].upper())) in tv:
   out.write(line); n+=1
   if n==1000: break
if n<1000: raise SystemExit(f'only {n} transversion markers available')
print(f'PASS: selected {n} deterministic Civán TV markers')
PY
python3 scripts/ecotype_pca_v2/10_call_ancient_fixed_markers.py --config "$RICE_PCA_CONFIG" --bam "$BAM" --fixed-snp "$FIXED" --sample LV7008416379 --panel C --library-type pooled_mixed --track TV --out-dir "$RICE_PCA_ATTEMPT_DIR/calls"
python3 scripts/ecotype_pca_v2/11_build_ancient_callability.py --config "$RICE_PCA_CONFIG" --fixed-snp "$FIXED" --calls "LV7008416379=$RICE_PCA_ATTEMPT_DIR/calls/LV7008416379.C.pooled_mixed.TV.calls.txt" --panel C --library-type pooled_mixed --track TV --out "$RICE_PCA_ATTEMPT_DIR/civan.callability.tsv"
printf 'metric\tvalue\nprototype\tCivan pooled_mixed TV\nmarker_n\t1000\nbam_interpretation\tpooled_mixed_capture_plus_shotgun\nslurm_job_id\t%s\n' "$SLURM_JOB_ID" > "$RICE_PCA_ATTEMPT_DIR/civan.prototype.summary.tsv"
echo "PASS: 50_civan_fixed_marker_prototype"
