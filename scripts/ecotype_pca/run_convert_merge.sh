#!/usr/bin/env bash
# ECOTYPE_PCA_PANEL.md 3.2节step2: 29M_3k(PLINK) -> EIGENSTRAT，再与6.7M_720求交集
# 需要: EIGENSOFT的convertf/mergeit在PATH里，python3+pysam(check_ref.py用)
set -euo pipefail

DB=/home/scratch/yinmt202607/db
PANEL_29M="$DB/29M_3k"
PANEL_720="$DB/6.7M_720"
OUT="$DB/merged_29M3k_6M7_720"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# a) 29M_3k: PLINK bed/bim.orig/fam -> EIGENSTRAT
#    convertf靠文件后缀自动识别PACKEDPED格式(snp文件必须.pedsnp/.map/.bim结尾，indiv文件必须
#    .pedind/.ped结尾)，没有inputformat这个参数。我们的原始文件名不符合(bim.orig/.fam)，
#    所以先建好后缀正确的软链接(不复制大文件)，再让convertf读
#    用裸数字染色体版 NB_final_snp.bim.orig，不要用改过染色体名(chr01)的版本
cd "$PANEL_29M"
ln -sf NB_final_snp.bim.orig NB_final_snp.rawchrom.bim
ln -sf NB_final_snp.fam NB_final_snp.rawchrom.pedind
convertf -p "$SCRIPT_DIR/par.PLINK.EIGENSTRAT"

# b) 核对convertf转换后的.snp，确认A1/A2方向没有被转错
#    （转换前29M_3k是A2=REF/A1=ALT，见ECOTYPE_PCA_PANEL.md 3.1节）
python3 "$SCRIPT_DIR/check_ref.py" "$PANEL_29M/NB_final_snp.snp" snp 200

# c) mergeit求交集：按染色体+物理位置匹配，自动纠正A1/A2顺序颠倒，
#    剔除strand-ambiguous(A/T、C/G)位点
mkdir -p "$OUT"
mergeit -p "$SCRIPT_DIR/par.MERGE"

echo "done: merged panel written to $OUT (merged.geno/.snp/.ind)"
