#!/bin/bash
set -euo pipefail

# Convert Civan et al. 2019's core SNP matrix (VCF) to EIGENSTRAT format,
# via an intermediate PLINK bed/bim/fam step.
#
# WHY TWO STEPS: convertf does NOT support VCF as an input format. Confirmed
# 2026-08-11 against the real CONVERTF/README on the server (not guessed) --
# convertf only supports 5 formats: ANCESTRYMAP, EIGENSTRAT, PED, PACKEDPED,
# PACKEDANCESTRYMAP. So plink2 does the VCF -> PACKEDPED (bed/bim/fam)
# conversion first, then convertf does PACKEDPED -> EIGENSTRAT, the same
# convertf step already used successfully for the 29M_3k panel.
#
# Usage: bash convert_civan_vcf.sh <db_dir> <script_dir>
#   db_dir:     directory containing sativa-rufipogon_SNPs.vcf.gz
#               (e.g. /home/scratch/yinmt202607/db/paper1)
#   script_dir: directory containing par.CIVAN.PLINK.EIGENSTRAT
#               (same as db_dir if you downloaded both files there together)

DB_DIR="${1:?usage: convert_civan_vcf.sh <db_dir> <script_dir>}"
SCRIPT_DIR="${2:?usage: convert_civan_vcf.sh <db_dir> <script_dir>}"

cd "$DB_DIR"

echo "[step 1] plink2: VCF -> PLINK bed/bim/fam"
echo "  --max-alleles 2 drops multiallelic sites (EIGENSTRAT/convertf are"
echo "  biallelic-only; check_ref.py's vcf mode already skips these when"
echo "  sampling, so this is consistent with how we already checked the data)"
echo "  --vcf-half-call missing treats any half-called genotype (e.g. 0/.)"
echo "  as missing rather than erroring out"
echo "  --chr-set 12 no-xy no-mt: rice has 12 autosomes, no X/Y/MT in this"
echo "  panel -- IF plink2 rejects this flag's syntax, run"
echo "  'plink2 --help 2>&1 | grep -A 8 -- \"--chr-set\"' and report back,"
echo "  the exact modifier keywords can differ slightly by plink2 build"
plink2 \
  --vcf sativa-rufipogon_SNPs.vcf.gz \
  --max-alleles 2 \
  --vcf-half-call missing \
  --chr-set 12 no-xy no-mt \
  --make-bed \
  --out civan_snp

echo "[step 2] symlink .fam -> .pedind"
echo "  convertf requires the indiv file to end in .pedind (or be a full"
echo "  .ped file) -- .fam is not an accepted suffix on its own. .bim IS"
echo "  already an accepted suffix for the snp file (see CONVERTF/README),"
echo "  so it does not need renaming, unlike the .fam file here."
ln -sf civan_snp.fam civan_snp.pedind

echo "[step 3] convertf: PACKEDPED -> EIGENSTRAT"
convertf -p "$SCRIPT_DIR/par.CIVAN.PLINK.EIGENSTRAT"

echo "done."
echo "output: civan_snp.eigenstratgeno / civan_snp.snp / civan_snp.ind"
echo "sanity check: compare SNP count in civan_snp.snp against the paper's"
echo "claimed 2,365,188 -- if wildly different, the VCF may not be complete"
echo "(see ECOTYPE_PCA_PANEL.md 5.6 -- total VCF line count was never"
echo "confirmed, only spot-checked)"
