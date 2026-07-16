#!/usr/bin/env bash
set -euo pipefail

# Download a small, reusable rice test panel for method dry-runs.
# Source convention checked on 2026-07-14 from Ensembl Plants species pages
# and their linked FTP indexes.

ROOT_DIR="${1:-genomes/rice_test_panel}"

mkdir -p "${ROOT_DIR}"/{nipponbare_irgsp10,azucena_rs1,ir64_osir64rs1,mh63_rs2,n22_osn22rs2}

download_pair() {
  local outdir="$1"
  local fasta_url="$2"
  local gff_url="$3"

  (
    cd "${outdir}"
    curl -fL --retry 3 -C - -O "${fasta_url}"
    curl -fL --retry 3 -C - -O "${gff_url}"
  )
}

# Standard japonica reference. Ensembl Plants release 63.
download_pair \
  "${ROOT_DIR}/nipponbare_irgsp10" \
  "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-63/fasta/oryza_sativa/dna/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa.gz" \
  "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-63/gff3/oryza_sativa/Oryza_sativa.IRGSP-1.0.63.gff3.gz"

# Practical upland-side proxy. Ensembl Plants release 62.
# GFF3 filename is inferred from the same release-62 naming pattern used by the
# other PSRefSeq cultivar FTP indexes linked from the species page.
download_pair \
  "${ROOT_DIR}/azucena_rs1" \
  "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/fasta/oryza_sativa_azucena/dna/Oryza_sativa_azucena.AzucenaRS1.dna.toplevel.fa.gz" \
  "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/gff3/oryza_sativa_azucena/Oryza_sativa_azucena.AzucenaRS1.62.gff3.gz"

# Practical lowland-side proxy. Ensembl Plants release 62.
download_pair \
  "${ROOT_DIR}/ir64_osir64rs1" \
  "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/fasta/oryza_sativa_ir64/dna/Oryza_sativa_ir64.OsIR64RS1.dna.toplevel.fa.gz" \
  "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/gff3/oryza_sativa_ir64/Oryza_sativa_ir64.OsIR64RS1.62.gff3.gz"

# Indica-side cultivated rice reference-like complement. Ensembl Plants release 62.
download_pair \
  "${ROOT_DIR}/mh63_rs2" \
  "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/fasta/oryza_sativa_mh63/dna/Oryza_sativa_mh63.MH63RS2.dna.toplevel.fa.gz" \
  "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/gff3/oryza_sativa_mh63/Oryza_sativa_mh63.MH63RS2.62.gff3.gz"

# Stress/aus-side line for robustness checks. Ensembl Plants release 62.
# GFF3 filename is inferred from the same release-62 naming pattern used by the
# verified IR64 and MH63 cultivar FTP indexes linked from the species page.
download_pair \
  "${ROOT_DIR}/n22_osn22rs2" \
  "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/fasta/oryza_sativa_n22/dna/Oryza_sativa_n22.OsN22RS2.dna.toplevel.fa.gz" \
  "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/gff3/oryza_sativa_n22/Oryza_sativa_n22.OsN22RS2.62.gff3.gz"

printf '\nDownloaded into: %s\n' "${ROOT_DIR}"
