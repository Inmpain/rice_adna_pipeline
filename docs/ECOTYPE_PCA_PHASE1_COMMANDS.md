# Ecotype PCA Phase 1: population labels, matrix filtering, sample-specific subset PCA, leave-one-out

This runbook starts from Phase 0's completed panel-overlap census
(`docs/ECOTYPE_PCA_PHASE0_COMMANDS.md`) and the whole-genus
`besthit_oryza.irgsp.bam` set (`gene/results/ecotype_pca/bam_irgsp/`).

**Decision recorded here, not repeated elsewhere**: the ORSC-narrowed
read set (`target_orsc`, Phase 0's step 3) turned out to throw away
70-90% of a sample's besthit-kept reads, which is too aggressive for
already-thin ancient samples -- LV7008416294's ORSC-narrowed IRGSP
mapping only reached mapq20_nondup=47. Phase 1 onward uses the
whole-genus `bam_irgsp/*.besthit_oryza.irgsp.bam` set instead
(`mapping_summary.tsv` in that directory shows mapped reads in the
1600-11800 range across all 16 samples, much healthier).

## 0. Cluster infrastructure note -- read before submitting anything

`/itp` is mounted on `node01-node04` but **not on `node06`** (confirmed
via `srun -w <node> mount | grep itp` per node; `node05` was down at
verification time and untested). Any `sbatch`/`srun` job whose path
resolves through an itp-backed symlink fails instantly and silently
(`FAILED`, exit `1:0`, 0 seconds, empty log) if scheduled onto node06.
**Every command below that uses `sbatch` includes `--exclude=node05,node06`
-- do not drop it.** For interactive `srun`, export once per shell
session instead:

```bash
export SBATCH_EXCLUDE=node05,node06
```

Reference data actively read by mapping/PCA jobs (anything under `db/`)
should live as real local directories, not itp symlinks -- itp is only
safe for data exclusively accessed from the login node.

## 1. Download the Phase 1 scripts

```bash
cd /home/scratch/yinmt202607/gene/scripts
BASE="https://api.github.com/repos/Inmpain/rice_adna_pipeline/contents/scripts/ecotype_pca"
REF="codex/ecotype-pca-panel"

for f in build_29m3k_population_labels.py build_720_population_labels.py \
         build_civan_population_labels.py filter_panel_by_label.py \
         build_sample_panel_subset.py simulate_leaveoneout_projection.py \
         run_sample_panel_pca.sh summarize_projection_distances.py \
         pseudo_haploid_call.py merge_ancient_into_panel.py; do
  curl -s "${BASE}/${f}?ref=${REF}" \
    | python3 -c "import sys,json,base64; print(base64.b64decode(json.load(sys.stdin)['content']).decode())" \
    > "$f"
  chmod +x "$f"
done
```

Using the GitHub Contents API (not `raw.githubusercontent.com`) deliberately
-- the raw-content CDN cached stale versions for several minutes after a
push during this work, even with a cache-busting query string; the API
endpoint did not.

## 2. Population labels: 29M_3k

Reuses `docs/references/3k_rice_genomes_project/rice_line_metadata_20141029.xlsx`
(this branch now carries its own copy, not just `main`).

```bash
mkdir -p /home/scratch/yinmt202607/db/29M_3k/references
curl -s "https://api.github.com/repos/Inmpain/rice_adna_pipeline/contents/docs/references/3k_rice_genomes_project/rice_line_metadata_20141029.xlsx?ref=codex/ecotype-pca-panel" \
  | python3 -c "import sys,json,base64; open('/home/scratch/yinmt202607/db/29M_3k/references/rice_line_metadata_20141029.xlsx','wb').write(base64.b64decode(json.load(sys.stdin)['content']))"

python3 build_29m3k_population_labels.py \
  --ind /home/scratch/yinmt202607/db/29M_3k/NB_final_snp.ind \
  --metadata-xlsx /home/scratch/yinmt202607/db/29M_3k/references/rice_line_metadata_20141029.xlsx \
  --out /home/scratch/yinmt202607/db/29M_3k/NB_final_snp.labeled.ind \
  --report /home/scratch/yinmt202607/db/29M_3k/NB_final_snp.label_report.tsv
```

Expect `3000/3024` matched; label distribution: IND 1743, TRJ 388, TEJ
319, AUS 215, INTERMEDIATE_TYPE 135, JAPONICA_UNSPEC 132, ARO 68, UNK 24.
`INTERMEDIATE_TYPE`/`JAPONICA_UNSPEC` are kept (not one of the six clean
codes, not folded into ADM without evidence) but naturally excluded from
`poplistname` since they're not in the standard label set.

## 3. Population labels: 6.7M_720

```bash
python3 build_720_population_labels.py \
  --ind /home/scratch/yinmt202607/db/6.7M_720/asn720.6m.ind \
  --pop-fam /home/scratch/yinmt202607/db/asn720data/asn720.pop.fam \
  --metadata-xlsx /home/scratch/yinmt202607/db/29M_3k/references/rice_line_metadata_20141029.xlsx \
  --out /home/scratch/yinmt202607/db/6.7M_720/asn720.6m.labeled.ind \
  --report /home/scratch/yinmt202607/db/6.7M_720/asn720.6m.label_report.tsv
```

Expect `718/720` matched (410 via direct `asn720.pop.fam` ID match, 308
of 310 `_merged`-suffix samples via the same 3K RGP metadata used in
step 2). `db/wild_rice_pangenome_README.txt` is empty on the server --
whether `OrA-OrF` matches the Nature 2025 wild-rice-pangenome paper's
`Or-Ia/Or-Ib/Or-II/Or-IIIa/Or-IIIb` naming is unresolved; labels are
passed through verbatim, not assumed equivalent.

## 4. Population labels: Civáň

Needs `Table_S1.csv` (sample metadata) and `Table_S2.csv` (bridges wild
samples' `W####` names to the `ERR......` SRA-run-accession IDs actually
used in `civan_snp.ind` -- confirmed 1:1 sequential, e.g. `W0101` ->
`ERR068593`), both already on the server under `db/paper1/`.

```bash
python3 build_civan_population_labels.py \
  --ind /home/scratch/yinmt202607/db/paper1/civan_snp.ind \
  --table-s1 /home/scratch/yinmt202607/db/paper1/Table_S1.csv \
  --table-s2 /home/scratch/yinmt202607/db/paper1/Table_S2.csv \
  --out /home/scratch/yinmt202607/db/paper1/civan_snp.labeled.ind \
  --report /home/scratch/yinmt202607/db/paper1/civan_snp.label_report.tsv
```

Expect `1055/1056` matched. Label distribution: indica 283, TRJ-style
"japonica (tropical)" 80 (written as `japonica_(tropical)`, whitespace
sanitized -- see step 4a), "japonica (temperate)" 51, aromatic 34,
unqualified japonica 23, O._rufipogon 456, four wild-outgroup singletons
(O._meridionalis/glaberrima/barthii/longistaminata, 1 each), UNK 1.

**4a. Known trap, already fixed in the script but worth knowing**: raw
`Group`/`Species` values from `Table_S1.csv` can contain literal spaces
("japonica (tropical)", "O. rufipogon") -- EIGENSTRAT `.ind` is
whitespace-delimited, so an unsanitized label breaks parsing for every
downstream consumer (this project's own `filter_panel_by_label.py`
hard-failed on exactly this on a real run, and smartpca itself would
likely mis-split it too). `build_civan_population_labels.py`'s
`sanitize_label()` already collapses internal whitespace to `_` --
nothing to do here, just don't reintroduce raw values if this script is
ever modified.

## 5. Drop unclassifiable (UNK) individuals from all three panels

`-lsqproject` projects any individual left in `.ind` but omitted from
`poplistname` -- leaving UNK samples in does not make them disappear
from the final plot, only physically removing them from both `.ind` and
`.eigenstratgeno`/`.geno` does. Each panel's genotype matrix is large
(29M_3k ~90GB), so run via `sbatch`, not interactively.

```bash
mkdir -p /home/scratch/yinmt202607/gene/results/ecotype_pca/loo_smoke

sbatch -p comp --exclude=node05,node06 -c 1 --mem 4G -t 12:00:00 \
  -J filter_29m3k -o /home/scratch/yinmt202607/db/29M_3k/filter_29m3k.%j.log \
  --wrap="python3 /home/scratch/yinmt202607/gene/scripts/filter_panel_by_label.py \
    --ind /home/scratch/yinmt202607/db/29M_3k/NB_final_snp.labeled.ind --drop-label UNK \
    --ind-out /home/scratch/yinmt202607/db/29M_3k/NB_final_snp.filtered.ind \
    --geno-in /home/scratch/yinmt202607/db/29M_3k/NB_final_snp.eigenstratgeno \
    --geno-out /home/scratch/yinmt202607/db/29M_3k/NB_final_snp.filtered.eigenstratgeno"

sbatch -p comp --exclude=node05,node06 -c 1 --mem 2G -t 02:00:00 \
  -J filter_720 -o /home/scratch/yinmt202607/db/6.7M_720/filter_720.%j.log \
  --wrap="python3 /home/scratch/yinmt202607/gene/scripts/filter_panel_by_label.py \
    --ind /home/scratch/yinmt202607/db/6.7M_720/asn720.6m.labeled.ind --drop-label UNK \
    --ind-out /home/scratch/yinmt202607/db/6.7M_720/asn720.6m.filtered.ind \
    --geno-in /home/scratch/yinmt202607/db/6.7M_720/asn720.6m.geno \
    --geno-out /home/scratch/yinmt202607/db/6.7M_720/asn720.6m.filtered.geno"

sbatch -p comp --exclude=node05,node06 -c 1 --mem 2G -t 02:00:00 \
  -J filter_civan -o /home/scratch/yinmt202607/db/paper1/filter_civan.%j.log \
  --wrap="python3 /home/scratch/yinmt202607/gene/scripts/filter_panel_by_label.py \
    --ind /home/scratch/yinmt202607/db/paper1/civan_snp.labeled.ind --drop-label UNK \
    --ind-out /home/scratch/yinmt202607/db/paper1/civan_snp.filtered.ind \
    --geno-in /home/scratch/yinmt202607/db/paper1/civan_snp.eigenstratgeno \
    --geno-out /home/scratch/yinmt202607/db/paper1/civan_snp.filtered.eigenstratgeno"

squeue -u "$USER"
```

Verify all three: `wc -l *.filtered.ind` should read 3000 / 718 / 1055,
and each job's log should end with `first-row width check passed`.
**All new scripts (step 6 onward) use the `.filtered.*` files -- never
the unfiltered originals.**

## 6. Smoke test: sample-specific subset PCA + leave-one-out, Civáň panel

Confirms the whole chain end-to-end and validates the REF/ALT 0/2
encoding direction before trusting `pseudo_haploid_call.py`'s output on
real data (see that script's docstring point 4 and
`docs/ECOTYPE_PCA_EXECUTION_PLAN.md` section 2.1/6).

```bash
cd /home/scratch/yinmt202607/gene/scripts
mkdir -p /home/scratch/yinmt202607/gene/results/ecotype_pca/loo_smoke
cd /home/scratch/yinmt202607/gene/results/ecotype_pca/loo_smoke

SAMPLE=LV7008416379   # best panel_overlap.tsv coverage of the 16 samples
HELD_OUT=$(awk '$3=="indica"{print $1; exit}' /home/scratch/yinmt202607/db/paper1/civan_snp.filtered.ind)

# ① + ②: pseudo-haploid calling (TV track, default)
python3 ../../../scripts/pseudo_haploid_call.py \
  --bam /home/scratch/yinmt202607/gene/results/ecotype_pca/bam_irgsp/${SAMPLE}.besthit_oryza.irgsp.bam \
  --panel-snp /home/scratch/yinmt202607/db/paper1/civan_snp.snp \
  --out ${SAMPLE}.civan.pseudohap.txt \
  --report ${SAMPLE}.civan.pseudohap.report.tsv
# Expect called=147 for LV7008416379 -- must match summarize_panel_overlap.py's
# callable_tv_high for this sample+panel exactly (independent cross-check).

# ②.5: shrink the panel to just this sample's 147 covered rows
python3 ../../../scripts/build_sample_panel_subset.py \
  --panel-snp /home/scratch/yinmt202607/db/paper1/civan_snp.snp \
  --panel-geno /home/scratch/yinmt202607/db/paper1/civan_snp.filtered.eigenstratgeno \
  --ancient-calls ${SAMPLE}.civan.pseudohap.txt \
  --out-panel-snp ${SAMPLE}.civan.subset.snp \
  --out-panel-geno ${SAMPLE}.civan.subset.eigenstratgeno \
  --out-ancient-calls ${SAMPLE}.civan.subset.calls.txt \
  --report ${SAMPLE}.civan.subset.report.tsv

# leave-one-out: mask $HELD_OUT (a known-indica sample) to LV7008416379's
# coverage pattern, re-call, then subset using --mask-from so both
# subsets align to the exact same 147 rows (the held-out individual has
# its own real missingness on top of the mask -- 24 of 147 sites, in
# this run -- naive per-file subsetting would misalign the two panels)
python3 ../../../scripts/simulate_leaveoneout_projection.py \
  --panel-geno /home/scratch/yinmt202607/db/paper1/civan_snp.filtered.eigenstratgeno \
  --panel-ind /home/scratch/yinmt202607/db/paper1/civan_snp.filtered.ind \
  --held-out-sample "$HELD_OUT" \
  --mask-from ${SAMPLE}.civan.pseudohap.txt \
  --seed 0 \
  --out ${HELD_OUT}.simulated.civan.pseudohap.txt \
  --report ${HELD_OUT}.simulated.civan.report.tsv

python3 ../../../scripts/build_sample_panel_subset.py \
  --panel-snp /home/scratch/yinmt202607/db/paper1/civan_snp.snp \
  --panel-geno /home/scratch/yinmt202607/db/paper1/civan_snp.filtered.eigenstratgeno \
  --ancient-calls ${HELD_OUT}.simulated.civan.pseudohap.txt \
  --mask-from ${SAMPLE}.civan.pseudohap.txt \
  --out-panel-snp ${HELD_OUT}.simulated.civan.subset.snp \
  --out-panel-geno ${HELD_OUT}.simulated.civan.subset.eigenstratgeno \
  --out-ancient-calls ${HELD_OUT}.simulated.civan.subset.calls.txt \
  --report ${HELD_OUT}.simulated.civan.subset.report.tsv

diff ${SAMPLE}.civan.subset.snp ${HELD_OUT}.simulated.civan.subset.snp \
  && echo "row-aligned, safe to merge"

# ③: relabel the held-out individual so it doesn't help build its own
# axis, then merge both the real ancient sample and the simulated
# individual into the same subset panel as two new columns
awk -v s="$HELD_OUT" '{if($1==s){$3="LOO_HELDOUT_EXCLUDED"}; print}' \
  /home/scratch/yinmt202607/db/paper1/civan_snp.filtered.ind > civan_snp.loo_test.ind

python3 ../../../scripts/merge_ancient_into_panel.py \
  --panel-geno ${SAMPLE}.civan.subset.eigenstratgeno \
  --panel-ind civan_snp.loo_test.ind \
  --calls ${SAMPLE}=${SAMPLE}.civan.subset.calls.txt \
          ${HELD_OUT}_LOOSIM=${HELD_OUT}.simulated.civan.subset.calls.txt \
  --ancient-poplabel Ancient \
  --out-geno merged.civan.subset.eigenstratgeno \
  --out-ind merged.civan.subset.ind

# ④: smartpca (found at ~/software/EIG/bin/smartpca on this cluster)
awk '{print $3}' merged.civan.subset.ind | sort -u \
  | grep -vE "^(Ancient|LOO_HELDOUT_EXCLUDED)$" > poplistname.civan_loo_test.txt

cat > par.civan_loo_test << EOF
genotypename:    $(pwd)/merged.civan.subset.eigenstratgeno
snpname:         $(pwd)/${SAMPLE}.civan.subset.snp
indivname:       $(pwd)/merged.civan.subset.ind
evecoutname:     $(pwd)/civan_loo_test.evec
evaloutname:     $(pwd)/civan_loo_test.eval
poplistname:     $(pwd)/poplistname.civan_loo_test.txt
lsqproject:      YES
numoutevec:      10
numoutlieriter:  0
numchrom:        12
numthreads:      2
EOF

smartpca -p par.civan_loo_test > smartpca.civan_loo_test.log 2>&1
```

**Expected result (2026-08-13 run, exact numbers)**: `${HELD_OUT}_LOOSIM`
(the simulated 147-site individual) projects to essentially the same
coordinates as `$HELD_OUT`'s own full-genotype row (labeled
`LOO_HELDOUT_EXCLUDED` in the same `.evec`) -- PC1 differed by 0.0006,
PC2 by 0.0001 in this run. This is the direct evidence that the 0/2
REF/ALT encoding is correct end to end; a reversed encoding would
project systematically toward the *opposite* population, not just
noisier. `LV7008416379` itself projected to PC1=0.0283, PC2=-0.0074 --
nearest real population by distance (excluding wild-outgroup singletons
with n<5, see `summarize_projection_distances.py --min-pop-size`):
aromatic (0.0083), then japonica variants, far from indica/aus/O._rufipogon.
This is one sample, 147 sites, a single pseudo-haploid draw -- a smoke
result, not a final call; the bootstrap uncertainty quantification
(execution plan section 6) hasn't been run yet.

```bash
python3 ../../../scripts/summarize_projection_distances.py \
  --evec civan_loo_test.evec --sample LV7008416379 \
  --out LV7008416379.civan.nearest.tsv
```

## 7. Scale out: 16 samples x 3 panels

`run_sample_panel_pca.sh` wraps steps ①-④ (not the leave-one-out
step -- that was a one-time pipeline validation, not part of the batch)
into one command per sample-panel pair, skipping any combination whose
`.evec` already exists.

**2026-08-15 update**: the civan block below now passes
`REFERENCE_LABELS_FILE` (9th arg,
`civan_domesticated_reference_labels.txt`) so all 16 civan combinations
get the wild-rice-excluded-from-axis fix (commit `cd5de6a`) from the
start, instead of needing a second corrective pass later. This fix was
verified on LV7008416379 alone before being applied here at scale (see
`docs/ECOTYPE_PCA_PANEL_QC_DESIGN.md` section 0) -- the qualitative
nearest-population ranking held, only the raw PC coordinates/distances
shifted (not directly comparable across the two axis-builder versions,
see that section for why). The 29m3k and 720 blocks are unchanged (no
`REFERENCE_LABELS_FILE`): 29m3k is 100% cultivated already so this
doesn't apply, and 720's correct reference set is still an open,
separate design question (section 3 of the QC design doc).

**Standing caveat, still true when this batch is run**: none of the
three panels have had MAF/LD pruning applied, and each ancient sample
still gets its own independently-subset marker set (the
"reference-first" redesign in `ECOTYPE_PCA_PANEL_QC_DESIGN.md` section
1 has not landed) -- per that doc's section 6 item 4, this batch's
results may need to be re-run once that redesign lands. Run this batch
as a first look, not a final one.

```bash
cd /home/scratch/yinmt202607/gene/scripts
mkdir -p /home/scratch/yinmt202607/gene/results/ecotype_pca/pca_runs

SAMPLES="LV6000619499 LV6000619917 LV6000620016 LV6000620032 LV6000620166 LV6000620172 LV6000654686 LV6000654698 LV7008416272 LV7008416280 LV7008416294 LV7008416329 LV7008416339 LV7008416349 LV7008416379 LV7008416407"
BAM_DIR=/home/scratch/yinmt202607/gene/results/ecotype_pca/bam_irgsp
OUT_DIR=/home/scratch/yinmt202607/gene/results/ecotype_pca/pca_runs

for sample in $SAMPLES; do
  BAM="${BAM_DIR}/${sample}.besthit_oryza.irgsp.bam"

  sbatch -p comp --exclude=node05,node06 -c 2 --mem 4G -t 02:00:00 \
    -J pca_civan -o "${OUT_DIR}/${sample}.civan.log" \
    --wrap="./run_sample_panel_pca.sh ${sample} civan ${BAM} \
      /home/scratch/yinmt202607/db/paper1/civan_snp.snp \
      /home/scratch/yinmt202607/db/paper1/civan_snp.filtered.eigenstratgeno \
      /home/scratch/yinmt202607/db/paper1/civan_snp.filtered.ind \
      ${OUT_DIR} TV \
      /home/scratch/yinmt202607/gene/scripts/civan_domesticated_reference_labels.txt"

  sbatch -p comp --exclude=node05,node06 -c 2 --mem 4G -t 02:00:00 \
    -J pca_29m3k -o "${OUT_DIR}/${sample}.29m3k.log" \
    --wrap="./run_sample_panel_pca.sh ${sample} 29m3k ${BAM} \
      /home/scratch/yinmt202607/db/29M_3k/NB_final_snp.snp \
      /home/scratch/yinmt202607/db/29M_3k/NB_final_snp.filtered.eigenstratgeno \
      /home/scratch/yinmt202607/db/29M_3k/NB_final_snp.filtered.ind \
      ${OUT_DIR} TV"

  sbatch -p comp --exclude=node05,node06 -c 2 --mem 4G -t 02:00:00 \
    -J pca_720 -o "${OUT_DIR}/${sample}.720.log" \
    --wrap="./run_sample_panel_pca.sh ${sample} 720 ${BAM} \
      /home/scratch/yinmt202607/db/6.7M_720/asn720.6m.snp \
      /home/scratch/yinmt202607/db/6.7M_720/asn720.6m.filtered.geno \
      /home/scratch/yinmt202607/db/6.7M_720/asn720.6m.filtered.ind \
      ${OUT_DIR} TV"
done

squeue -u "$USER"
```

After all 48 finish, summarize every one:

```bash
cd "$OUT_DIR"
for evec in *.evec; do
  prefix="${evec%.evec}"
  sample="${prefix%%.*}"
  python3 /home/scratch/yinmt202607/gene/scripts/summarize_projection_distances.py \
    --evec "$evec" --sample "$sample" --out "${prefix}.nearest.tsv"
done

head -1 $(ls *.nearest.tsv | head -1) > all_samples.nearest.tsv
for f in *.nearest.tsv; do [[ "$f" == all_samples* ]] && continue; tail -n +2 "$f"; done >> all_samples.nearest.tsv
column -t -s $'\t' all_samples.nearest.tsv | less -S
```

Status as of this doc's last update (2026-08-15): not yet submitted.
Before submitting, `docs/ECOTYPE_PCA_PANEL_QC_DESIGN.md` section 0's
masked-LOO classification-power validation (32 individuals, aromatic/
aus/indica/japonica, on LV7008416379's own 147-site Civáň mask) is worth
reading first -- aromatic/aus/indica each classified with 100% recall
and aromatic had 88.9% precision (8/9 predicted-aromatic calls were
truly aromatic, 1 was a mislabeled japonica), but the plain "japonica" label
itself only had 37.5% recall (mostly confused with
japonica_(tropical), not with aromatic) -- so this marker set
discriminates aromatic well but has weak resolution within the
japonica-family subdivisions. Doesn't block running this batch, just
context for how much weight to put on any single sample's civan
"closest to X" result once this batch completes.
