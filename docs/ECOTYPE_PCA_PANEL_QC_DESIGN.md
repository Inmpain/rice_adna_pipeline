# Ecotype PCA panel QC design: MAF, LD pruning, reference-first architecture

**Status as of 2026-08-13**: this is a design/decision doc, not a
completed runbook. Only item 0 below has been implemented and pushed.
Everything else is a recorded plan -- if you are picking this up cold,
read this whole file before touching any panel's genotype matrix, then
re-read `docs/ECOTYPE_PCA_PANEL.md` and `docs/ECOTYPE_PCA_EXECUTION_PLAN.md`
for the broader pipeline this feeds into.

## Why this doc exists

The three panels (`29M_3k`, `6.7M_720`, Civáň) have so far only had
one filter applied: dropping fully-unclassifiable (`UNK`) samples
(`docs/ECOTYPE_PCA_PHASE1_COMMANDS.md` step 5). **None of the three
have had MAF filtering or LD pruning applied for PCA purposes.** The
user asked (2026-08-13) for a plain description of the three panels to
send to GPT for external review specifically on this gap. GPT's review
raised real, concrete design problems -- this doc records what GPT
said, what I independently verified against our own repo/code (not
just accepted), what's already fixed, and what's still open, so none
of it gets lost to context compaction.

## 0. Fixed already: Civáň panel PCA axes were being built partly from wild rice

**Commit `cd5de6a`, pushed to `codex/ecotype-pca-panel`.**

`run_sample_panel_pca.sh` used to build smartpca's `poplistname` from
*every* non-`Ancient` label present in the merged `.ind` file. On the
Civáň panel that silently included the 456 `O._rufipogon` wild samples
plus 4 singleton wild outgroups (`O._meridionalis`/`glaberrima`/
`barthii`/`longistaminata`) as axis-building references, not just
projected individuals. This is backwards from Civáň et al. 2019's own
method: build axes from the 595 domesticated accessions only, project
wild rice onto them (the paper's stated reason is wild rice's much
higher missingness -- letting it into axis-building would let coverage
artifacts, not just ancestry, shape the PCs).

**Fix**: `run_sample_panel_pca.sh` now takes an optional 9th argument,
`REFERENCE_LABELS_FILE` -- when given, `poplistname` is built from
exactly those labels (validated: the script hard-fails if any listed
label has zero individuals in the merged `.ind`, e.g. a typo). Omit the
argument to keep the old include-everything behavior, which is correct
for `29M_3k` (100% cultivated already, no wild-vs-domesticated split to
worry about) but is a **separate, not-yet-decided** design question for
`6.7M_720` (see section 3).

New file `scripts/ecotype_pca/civan_domesticated_reference_labels.txt`,
the 6 domesticated labels for the Civáň panel specifically:
```
indica
aus
aromatic
japonica
japonica_(tropical)
japonica_(temperate)
```

**Not yet done**: re-running LV7008416379's Civáň projection with this
fix to see whether "closest to aromatic" survives. Commands are already
given to the user (see chat), result not yet returned as of this doc's
last edit -- check `docs/ECOTYPE_PCA_PANEL.md`'s 📍 handoff block or ask
before re-deriving.

## 1. Core architectural principle GPT's review established (not yet implemented)

**Reference-first, frozen marker set.** For each panel, independently:

1. Fix which individuals count as "reference" for that panel (all of
   it for `29M_3k`; domesticated-only for Civáň, per section 0; TBD for
   `6.7M_720`, section 3).
2. Do coordinate/allele harmonization, biallelic-SNP check, QC --
   computed **only on the reference individuals**, never on ancient
   samples and never on non-reference modern individuals (e.g. Civáň's
   wild rice).
3. MAF filter -- again, allele frequency computed only within the
   reference set.
4. LD prune the reference set to a physical-window target.
5. **Freeze** this pruned SNP list as the panel's one-and-only PCA
   marker set (e.g. `panelX.pca.prune.in`).
6. Run smartpca **once per panel**, `poplistname` = the reference
   labels, covering **all ancient samples in a single run** (each as
   its own column, with 9/missing at whatever fraction of the frozen
   marker set that sample's reads didn't cover) plus any non-reference
   modern individuals (Civáň's wild rice) also just projected.
7. Only after that: masked-modern leave-one-out validation, on the
   same frozen marker set.

**Why this replaces the current per-sample-subset design**:
`build_sample_panel_subset.py` (used by `run_sample_panel_pca.sh` step
2) currently shrinks the panel to *each ancient sample's own covered
SNPs* before running smartpca. `-lsqproject` doesn't need this --
it's specifically designed to project a sparse/missing individual
against a fixed full marker set, treating uncovered positions as
missing. Per-sample subsetting means two different ancient samples get
two different marker sets, so their "PC1"/"PC2" are technically
different coordinate systems even though both get called PC1/PC2 in
the `.evec` output -- comparing them directly (e.g. building a temporal
trajectory across samples) is not currently valid. **This is the
single biggest pending redesign** and touches `build_sample_panel_subset.py`,
`merge_ancient_into_panel.py`'s call pattern, and `run_sample_panel_pca.sh`'s
per-sample-loop structure in section 7 of `ECOTYPE_PCA_PHASE1_COMMANDS.md`.
Not started. The existing per-sample-subset smoke-test results
(LV7008416379 -> aromatic) are still useful as a first-pass sanity
check, just not as a final, cross-sample-comparable result.

## 2. Panel 1: 29M_3k (cultivated, PCA-A)

- Citation, confirmed against `docs/LITERATURE.md` section 3 (main
  branch, not re-derived from GPT): Wang et al. 2018, *Nature*
  557:43-49, "Genomic variation in 3,010 diverse accessions of Asian
  cultivated rice", DOI 10.1038/s41586-018-0063-9; plus two GigaScience
  2014 papers. 3024 samples, ~29.6M biallelic SNP, PLINK format,
  Nipponbare MSU7/IRGSP1.0 coordinates.
- No wild-vs-domesticated axis question -- IND/AUS/ARO/TRJ/TEJ/ADM are
  all peer-level cultivated groups, all can act as axis-building
  reference (no exclusion list needed here, unlike Civáň).
- **Still needs its own independent MAF + LD pruning pass** -- this
  hasn't been touched at all beyond UNK removal. Not started.
- GPT's suggested primary parameters (see section 4 for the caveats on
  treating these as defaults rather than final answers):
  `--maf 0.01` computed within the 3000-sample cultivated reference,
  `--indep-pairwise 100kb 10 0.2` as primary LD-pruning window.

## 3. Panel 2: 6.7M_720 (wild-dominant, PCA-B) -- most urgent, least understood

- Citation, **independently verified via WebSearch (not just accepted
  from GPT)**: Wang H, Vieira FG, Crawford JE, Chu C, Nielsen R (2017).
  "Asian wild rice is a hybrid swarm with extensive gene flow and
  feralization from domesticated rice." *Genome Research* 27(6):1029-1038.
  DOI: 10.1101/gr.204800.116. PMID: 28385712. Confirmed via two
  independent WebSearch queries: title/authors/journal match, and the
  paper's own sample composition (203 domesticated + 435 wild = 638
  total) matches what GPT stated.
- **What we were NOT able to independently verify**: GPT's specific
  methods claim that the original paper computed allele frequencies via
  `ANGSD -doMaf` and reduced to markers by randomly keeping one variable
  site per 5kb window, arriving at 60,722 markers for the population-
  structure/PCA analysis. This is plausible and consistent with the
  638-sample count, but we have not pulled the paper's Methods section
  ourselves to confirm the exact numbers -- treat as GPT-relayed,
  not confirmed, until someone reads the actual paper text.
- **The data we actually have is NOT that 60,722-marker set.** Per
  `docs/ECOTYPE_PCA_PANEL.md` section 1.2, our `asn720.6m.*` files are
  720 samples (not 638) and ~6.7M SNP sites (not 60,722) -- the
  panel-overlap docs describe this as "the first author later increased
  the site density and sent it to us directly," i.e. **a different,
  custom re-processing of (probably, not confirmed) the same underlying
  sequencing data, not the paper's own published analysis-ready
  matrix.** Its MAF/missingness/LD-pruning status is therefore
  genuinely unknown to us, not just "not yet filtered by us."
- **Format correction to GPT's suggested commands**: GPT's audit script
  proposal assumes `plink2 --pfile rice720 ...` (PLINK2 native
  `.pgen/.pvar/.psam`). Our `6.7M_720` panel is **EIGENSTRAT**
  (`.geno`/`.ind`/`.snp` text files, confirmed in
  `docs/ECOTYPE_PCA_PANEL.md` 1.2), not a PLINK pfile. Before any of
  GPT's PLINK2 commands can run as-written, need either:
  - `convertf` (EIGENSOFT) to convert EIGENSTRAT -> PACKEDANCESTRYMAP,
    or
  - `plink2 --eigfile PREFIX ...` / separate `--eiggeno`/`--eigind`/
    `--eigsnp` flags (PLINK2 does support importing EIGENSOFT binary
    directly -- but check whether our `.geno` is plain-text EIGENSTRAT
    or already packed binary before assuming which import path
    applies; PHASE1 commands never established this, it needs a fresh
    `file`/`head -c 100` check on the real files).
  Not started -- confirm PLINK2 is even installed on the cluster
  (`module avail plink2` / `which plink2`) before writing final
  commands.
- **Planned audit** (not yet run), in this order:
  1. `file`/`head -c 100` check on `asn720.6m.geno` to confirm text vs
     binary EIGENSTRAT, so the right convertf/plink2 import path can be
     written.
  2. Allele-frequency (MAF) distribution across the whole 6.7M sites,
     binned: `0`, `<0.001`, `0.001-0.005`, `0.005-0.01`, `0.01-0.05`,
     `0.05-0.10`, `>0.10`. With 720 diploid individuals, singleton MAF
     is theoretically `1/1440 ≈ 0.000694`, doubleton `2/1440 ≈ 0.00139`
     -- if the distribution shows a hard floor near one of these values
     rather than a smooth tail, that's evidence some MAF filter was
     already applied upstream by the author.
  3. Per-site and per-individual missingness distributions -- flagged
     as especially important here since the panel mixes ID styles
     (`ERR0685xx`-style ENA run accessions vs `B0xx_merged`-style 3K-RGP
     names, per `docs/ECOTYPE_PCA_PANEL.md` 1.2), suggesting the 720
     samples were pooled from different sequencing batches with
     potentially very different depth/missingness -- if e.g. one `OrX`
     population label has systematically higher missingness than
     others, that's a real confound for any PCA built on this panel.
  4. Marker physical-spacing distribution (adjacent-SNP distance per
     chromosome) -- if the original paper's "1 SNP per 5kb" thinning
     survived into this file, spacing would cluster near 5kb; 6.7M
     sites across a ~380Mb genome makes that look unlikely on its face
     (implies ~57 sites/kb genome-wide) but should be checked directly,
     not assumed.
  5. LD decay by physical-distance bin (pairwise r², binned `0-1kb`,
     `1-2kb`, `2-5kb`, `5-10kb`, `10-20kb`, `20-50kb`, `50-100kb`,
     `100-200kb`, `200-500kb`), on a handful of chromosomes (suggested:
     chr1, chr3, chr6, chr12) -- this determines what LD-pruning window
     is actually appropriate here, rather than assuming the `100kb`
     default that might make sense for `29M_3k`.
  6. Only after 1-5: decide between two candidate pruning routes --
     **Route A** (mimic the original paper's approach): keep at most 1
     variable site per 5kb window. **Route B** (standard PLINK LD
     pruning): `--maf 0.01 --geno 0.05 --indep-pairwise <window from
     step 5> 0.2`, with `r² < 0.1/0.2/0.5` and multiple window sizes as
     sensitivity checks, mirroring the `29M_3k` sensitivity-analysis
     approach in section 4.
- **Open, unresolved: which individuals should count as "reference" for
  axis-building on this panel.** Unlike Civáň (clean cultivated-vs-wild
  split, domesticated defines axes per the paper's own method), 720 is
  wild-dominant and the `OrA`-`OrF` labels' correspondence to any
  published ecotype/lineage scheme is still unconfirmed (per
  `docs/ECOTYPE_PCA_PANEL.md` 1.2's "OrA-OrF不是最终结论" caveat, itself
  still open). GPT suggested "selected modern wild reps" as axis
  reference but did not give a concrete selection rule -- do not
  implement `civan`-style domesticated-only exclusion here without
  first resolving what "reference" even means for this panel. Also:
  GPT recommended treating this panel's scientific scope as "cultivated-
  like vs wild-rice-related / wild population affinity" only, **not**
  giving ancient samples' `OrA`-`OrF` projections an ecotype
  interpretation until the label correspondence to Wang et al. 2017's
  own population definitions is confirmed -- carry this caveat forward
  into any results written up from this panel.

## 4. Panel 3: Civáň (bridging panel, PCA-C)

- Citation: confirmed independently already in `docs/LITERATURE.md`
  section 2.2 (main branch) before this GPT conversation happened --
  Civáň P, Ali S, Batista-Navarro R, Drosou K, Ihejieto C, Chakraborty
  D, Ray A, Gladieux P, Brown TA (2019). *Genome Biology and Evolution*
  11(3):832-843. DOI: 10.1093/gbe/evz039. PMID: 30793171. 1056 samples
  (595 cultivated: 283 indica/154 japonica/124 aus/34 aromatic + 461
  wild, mostly O. rufipogon), 2,365,188 biallelic SNP, IRGSP-1.0.
- Axis-building fix: done, see section 0.
- **Still needs MAF + LD pruning** -- the smoke test so far
  (`docs/ECOTYPE_PCA_PHASE1_COMMANDS.md` section 6) ran on the raw
  2.365M-SNP `.filtered.*` matrix (UNK removed, nothing else). Not
  started.
- GPT's caveat, worth keeping verbatim: the paper's own PCA used two
  *different* SNP sets for two different purposes -- an LD-pruned
  ~404k "CoreSNP" set to pick "typical" indica/japonica/aus/aromatic
  representative individuals, versus the full ~2.365M merged matrix
  (which the paper itself describes as not MAF/LD-clean -- lots of
  low-frequency variation in the cultivated side, high missingness in
  the wild subset) for the actual ancestry-painting analysis. **We have
  not independently confirmed these two specific numbers (404k core
  SNP count, and the exact QC steps behind it) against the paper's own
  text** -- flagged as GPT-relayed pending independent verification,
  same caveat as the 720-panel methods claim above.
- **Sensitivity analyses GPT recommended, not yet run**:
  - MAF: primary `0.01` (computed on the 595-sample domesticated
    reference only, not the full 1056 including wild), sensitivity
    `0.05`. Reasoning given: aromatic is already a small reference
    group (34 samples) -- an aggressive MAF cutoff risks losing exactly
    the low-frequency variation that carries information about small
    subgroups.
  - LD pruning: primary `--indep-pairwise 100kb 10 0.2`, sensitivity
    `50kb`/`100kb`/`200kb` at the same step/r² -- goal is checking
    whether "closest to aromatic" is stable across these choices, not
    picking one "correct" value.
  - **PC dimensionality**: don't stop at PC1-PC2. Paper apparently
    shows aromatic's separation from japonica is clearer once PC3 is
    included (2D PC1-PC2 alone partially conflates them) -- GPT
    recommended checking nearest-population ranking at PC1-3, PC1-5,
    and PC1-10, not just the PC1-2 pair `summarize_projection_distances.py`
    currently defaults to (`--num-pcs` already supports this, just
    hasn't been run with a value above 2).
  - **LOO validation upgrade**: current smoke test masks exactly one
    known individual and checks it lands near its own centroid. GPT's
    stronger suggestion: mask many individuals across each reference
    population (not just one), and report per-population accuracy /
    false-positive rate (e.g. "of N known-aromatic individuals masked
    to ancient-sample coverage, what fraction still rank aromatic as
    nearest, and what fraction of known-indica/japonica individuals
    incorrectly rank aromatic as nearest") -- this is a much stronger
    piece of evidence than a single masked individual's centroid
    distance, but requires running `simulate_leaveoneout_projection.py`
    across many held-out samples instead of one, which
    `simulate_leaveoneout_projection.py` already supports per-invocation
    (one `--held-out-sample` at a time) but has never been looped over
    a full population. Not started.

## 5. Standing caveat carried over, unrelated to MAF/LD but raised in the same GPT conversation: `ARO` label ≠ guaranteed fragrance

Civáň's own data: only 8 of their 34 aromatic accessions carry the
classic `badh2.1` fragrance-deletion allele; across the full 3K RGP set
of 76 "aromatic"-labeled accessions, only 26 carry it. **A sample
projecting nearest to the `aromatic`/`ARO` population label supports
"genome-wide affinity to the aromatic/circum-Basmati-related genetic
group," not "this rice was fragrant."** If/when an ancient sample's
PCA result gets written up as aromatic-affiliated, the fragrance claim
specifically needs an independent `BADH2` check, not an inference from
the population label alone. (This appears to be independently
GPT-relayed, not yet cross-checked against the paper's own text the way
the citation itself was -- same "verify before quoting as fact"
caveat as sections 3-4 above.)

## 6. Open TODOs, in the order they'd naturally get picked up

1. Get the LV7008416379-on-Civáň-domesticated-only-axes re-run result
   back (commands already given to the user, section 0) and compare
   against the old wild-inclusive result.
2. `6.7M_720` audit (section 3, steps 1-5) -- write the audit script,
   run it, decide Route A vs Route B pruning. This was explicitly
   flagged as the most urgent because its provenance is the least
   understood of the three panels.
3. MAF + LD pruning for `29M_3k` and Civáň (sections 2 and 4), using
   whatever window/MAF sensitivity grid gets decided.
4. The big one: reference-first architecture redesign (section 1) --
   rewrite `build_sample_panel_subset.py`/`merge_ancient_into_panel.py`/
   `run_sample_panel_pca.sh` so each panel gets ONE frozen pruned marker
   set and ONE smartpca run covering all 16 ancient samples together,
   replacing the current per-sample-subset-per-run design. This
   invalidates re-running section 7 of `ECOTYPE_PCA_PHASE1_COMMANDS.md`
   as currently written -- don't scale that out further until this
   redesign lands, or the 16-sample x 3-panel batch will need to be
   thrown away and re-run anyway.
5. Per-population LOO accuracy/false-positive-rate validation (section
   4's LOO upgrade), once the frozen marker set exists to run it on.
