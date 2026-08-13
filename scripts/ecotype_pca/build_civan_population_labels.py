#!/usr/bin/env python3
"""
Attach real population labels to civan_snp.ind from Civáň et al. 2019's
own Table_S1.csv sample metadata.

The panel's .ind file currently carries a single placeholder value,
"???", in its population column for all 1,056 samples -- this is
docs/ECOTYPE_PCA_PANEL.md section 3.2 todo item 2c.

ID FORMAT: civan_snp.ind mixes samples from sub-collections that entered
the source VCF with and without a pre-existing FID. Samples that had none
got plink2's default FID_IID self-duplication around a literal underscore
-- "B006" becomes "B006_B006", "IRIS_313-9986" (itself containing an
underscore) becomes "IRIS_313-9986_IRIS_313-9986" (see
docs/ECOTYPE_PCA_PANEL.md 1.3/5.6). Samples that already had a real FID
(confirmed on the real file: "ERR068593"-style IDs) were NOT doubled and
appear in civan_snp.ind exactly as in Table_S1.csv's Accession column.
Recovering the doubled form CANNOT be done by splitting on the first or
last underscore (that breaks for accessions like "IRIS_313-9986" which
contain underscores themselves) -- it must find the exact midpoint: the
string has odd length 2n+1, character n is '_', and the two n-character
halves are identical. recover_accession() below checks for that pattern
first and falls back to the sample ID unchanged (not doubled) when it
doesn't match, rather than hard-failing -- both forms are confirmed-real
cases here, not a guess.

LABEL SELECTION FOR WILD SAMPLES: Table_S1.csv's Group column is "-" for
every wild accession (confirmed on the real file) -- using that verbatim
would merge all wild samples under one meaningless "-" label regardless
of species, which is wrong, not just uninformative: of the 460 wild
samples resolved via the Table_S2 bridge below, 456 are O. rufipogon but
a handful are entirely different species (O. meridionalis, O.
glaberrima, O. barthii, O. longistaminata -- confirmed on the real
file). When Group is empty or "-", the label falls back to Species
instead, so these don't get silently pooled with O. rufipogon into the
same population for axis-building.

WILD-SAMPLE ID BRIDGE: the 461 wild samples are "plain"-form ERR-style
run accessions in civan_snp.ind, but Table_S1.csv identifies wild
accessions by the paper's own "W####" names (confirmed on the real
files: a first run against Table_S1.csv alone matched all 595
cultivated samples and exactly 0 of the 461 wild ones). Table_S2.csv
(nominally "chloroplast assembly QC") turns out to carry a direct,
1:1, sequential Accession<->"SRA dataset used" mapping for these same
W#### names (W0101->ERR068593, W0102->ERR068594, ... confirmed on the
real file) -- --table-s2 is optional specifically so this bridge can be
skipped if Table_S2.csv's coverage turns out incomplete (it has 1,645
rows against 1,825 total chloroplast genomes, a known partial-coverage
gap per docs/ECOTYPE_PCA_PANEL.md 1.3, so it may not cover all 461
nuclear-SNP wild samples either -- this run's [labels] summary reports
exactly how many it resolves).

WHY THIS SCRIPT DOES NOT STANDARDIZE THE LABEL, UNLIKE THE OTHER TWO
build_*_population_labels.py SCRIPTS: unlike NB_final_snp.ind's 8-value
xlsx column (fully scanned before that mapping was written), Civáň's
Table_S1.csv "Group" column values weren't known ahead of time because
its messy header (quoted fields containing literal newlines, 1,024
columns of which 1,013 are blank -- an Excel export artifact) initially
blocked reading it at all. The raw "Group" value is therefore used AS
THE LABEL VERBATIM (a real run now shows it's already clean: indica/aus/
aromatic/"japonica (tropical)"/"japonica (temperate)"/unqualified
"japonica", matching the paper's stated 283/124/34/80/51/23 split
exactly) and "Species" is carried in the report for context. A
standardization mapping to IND/AUS/ARO/TRJ/TEJ (folding the unqualified
"japonica" into JAPONICA_UNSPEC, matching NB_final_snp.ind's convention)
can be layered on as a followup now that the value list is confirmed.

Table_S1.csv was also feared (docs/ECOTYPE_PCA_PANEL.md 1.3, based on a
raw `wc -l`) to have 1,063 data rows against a paper-stated 1,056 -- a
proper CSV-aware read (this script) finds only 1,057, and the 6-line gap
was `wc -l` counting literal newlines embedded inside quoted header
cells as extra "lines", not a real data discrepancy. The single
remaining extra row is tolerated via last-row-wins on duplicate
Accession (reported, not hard-failed on, since this is an already-messy
file where a hard stop would block ever seeing the real numbers).

Usage:
  python3 build_civan_population_labels.py \\
    --ind /home/scratch/yinmt202607/db/paper1/civan_snp.ind \\
    --table-s1 /home/scratch/yinmt202607/db/paper1/Table_S1.csv \\
    --table-s2 /home/scratch/yinmt202607/db/paper1/Table_S2.csv \\
    --out /home/scratch/yinmt202607/db/paper1/civan_snp.labeled.ind \\
    --report /home/scratch/yinmt202607/db/paper1/civan_snp.label_report.tsv
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
import tempfile
from collections import Counter
from pathlib import Path


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--ind", required=True, help="path to civan_snp.ind")
    parser.add_argument("--table-s1", required=True, help="path to Civáň et al. 2019 Table_S1.csv")
    parser.add_argument(
        "--table-s2",
        default=None,
        help="optional path to Table_S2.csv, bridges wild-sample ERR-style IDs to Table_S1's W#### accessions",
    )
    parser.add_argument("--out", required=True, help="path for the new .ind with real labels")
    parser.add_argument("--report", required=True, help="path for a per-sample coverage TSV")
    parser.add_argument(
        "--unmapped-label",
        default="UNK",
        help="population label for samples with no match in Table_S1.csv (default: UNK)",
    )
    return parser.parse_args(argv)


def validate_args(args: argparse.Namespace) -> None:
    if not Path(args.ind).is_file():
        raise FileNotFoundError(f".ind file not found: {args.ind}")
    if not Path(args.table_s1).is_file():
        raise FileNotFoundError(f"Table_S1.csv not found: {args.table_s1}")
    if args.table_s2 is not None and not Path(args.table_s2).is_file():
        raise FileNotFoundError(f"Table_S2.csv not found: {args.table_s2}")


def sanitize_label(value: str) -> str:
    """Collapse internal whitespace to '_' so the label is one token.

    EIGENSTRAT .ind is whitespace-delimited (sample, sex, label) -- any
    consumer that splits on whitespace (smartpca itself, and this
    project's own filter_panel_by_label.py) will misparse a label
    containing a literal space. Confirmed on the real file: raw Group
    values like "japonica (tropical)" and Species values like
    "O. rufipogon" both contain one, and this was NOT caught until
    filter_panel_by_label.py hard-failed on a real run against the
    server's civan_snp.labeled.ind.
    """
    return re.sub(r"\s+", "_", value.strip())


def recover_accession(sample_id: str) -> tuple[str, str]:
    """Return (accession, method), method is 'doubled' or 'plain'."""
    n = len(sample_id)
    if n % 2 == 1:
        mid = n // 2
        if sample_id[mid] == "_" and sample_id[:mid] == sample_id[mid + 1 :]:
            return sample_id[:mid], "doubled"
    return sample_id, "plain"


def load_table_s1(path: Path) -> tuple[dict[str, tuple[str, str]], int, int]:
    """Return accession -> (species, group), total data rows, duplicate accession count."""
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        if len(header) < 3:
            raise ValueError(f"Table_S1.csv header has only {len(header)} columns, expected at least 3 (Accession, Species, Group)")
        if header[0].strip() != "Accession":
            raise ValueError(f"Table_S1.csv column 0 is {header[0]!r}, expected 'Accession'")
        if header[2].strip() != "Group":
            raise ValueError(f"Table_S1.csv column 2 is {header[2]!r}, expected 'Group'")

        lookup: dict[str, tuple[str, str]] = {}
        total_rows = 0
        duplicate_count = 0
        for row in reader:
            if not any(cell.strip() for cell in row):
                continue  # fully blank row
            total_rows += 1
            accession = row[0].strip()
            species = row[1].strip() if len(row) > 1 else ""
            group = row[2].strip() if len(row) > 2 else ""
            if not accession:
                continue
            if accession in lookup:
                duplicate_count += 1
            lookup[accession] = (species, group)
    return lookup, total_rows, duplicate_count


def load_table_s2(path: Path) -> dict[str, str]:
    """Return SRA-run-accession ("SRA dataset used" column) -> Table_S1 Accession."""
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        if len(header) < 3:
            raise ValueError(f"Table_S2.csv header has only {len(header)} columns, expected at least 3 (Accession, Species/group, SRA dataset used)")
        if header[0].strip() != "Accession":
            raise ValueError(f"Table_S2.csv column 0 is {header[0]!r}, expected 'Accession'")
        if header[2].strip() != "SRA dataset used":
            raise ValueError(f"Table_S2.csv column 2 is {header[2]!r}, expected 'SRA dataset used'")

        bridge: dict[str, str] = {}
        for row in reader:
            if not any(cell.strip() for cell in row):
                continue
            accession = row[0].strip()
            sra_id = row[2].strip() if len(row) > 2 else ""
            if not accession or not sra_id:
                continue
            # last-row-wins, matching load_table_s1's tolerance for this
            # already-messy file family rather than hard-failing
            bridge[sra_id] = accession
    return bridge


def load_ind(ind_path: Path) -> list[tuple[str, str, str]]:
    rows: list[tuple[str, str, str]] = []
    seen: set[str] = set()
    with ind_path.open() as handle:
        for line_no, line in enumerate(handle, start=1):
            fields = line.split()
            if len(fields) != 3:
                raise ValueError(
                    f"{ind_path}:{line_no}: expected 3 whitespace-separated "
                    f"fields (sample sex label), got {len(fields)}: {line!r}"
                )
            sample_id, sex, old_label = fields
            if sample_id in seen:
                raise ValueError(f"duplicate sample ID in {ind_path}: {sample_id!r}")
            seen.add(sample_id)
            rows.append((sample_id, sex, old_label))
    return rows


def write_atomic(path: Path, write_fn) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    os.close(fd)
    tmp_path = Path(tmp_name)
    try:
        with tmp_path.open("w", newline="") as handle:
            write_fn(handle)
        os.replace(tmp_path, path)
    except BaseException:
        tmp_path.unlink(missing_ok=True)
        raise


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        validate_args(args)
        table_s1, total_rows, dup_count = load_table_s1(Path(args.table_s1))
        table_s2 = load_table_s2(Path(args.table_s2)) if args.table_s2 else {}
        ind_rows = load_ind(Path(args.ind))

        report_rows = []
        new_ind_lines = []
        label_counts: Counter[str] = Counter()
        species_group_counts: Counter[tuple[str, str]] = Counter()
        method_counts: Counter[str] = Counter()
        bridge_used_count = 0
        n_matched = 0

        for sample_id, sex, _old_label in ind_rows:
            accession, method = recover_accession(sample_id)
            method_counts[method] += 1
            meta = table_s1.get(accession)
            bridged = False
            if meta is None and accession in table_s2:
                bridged_accession = table_s2[accession]
                meta = table_s1.get(bridged_accession)
                if meta is not None:
                    accession = bridged_accession
                    bridged = True
                    bridge_used_count += 1
            if meta is None:
                final_label = args.unmapped_label
                species, group = "", ""
            else:
                species, group = meta
                if group and group != "-":
                    final_label = sanitize_label(group)
                elif species:
                    # Table_S1.csv leaves Group as "-" for every wild
                    # accession (confirmed on the real file), which is not
                    # one label -- 456 of these rows are O. rufipogon but a
                    # handful are entirely different wild species
                    # (O. meridionalis / O. glaberrima / O. barthii /
                    # O. longistaminata, confirmed on the real file). Fall
                    # back to Species so these don't get silently merged
                    # into the same population label as O. rufipogon.
                    final_label = sanitize_label(species)
                else:
                    final_label = args.unmapped_label
                n_matched += 1
                species_group_counts[(species, group)] += 1
            label_counts[final_label] += 1
            report_rows.append(
                {
                    "sample_id": sample_id,
                    "accession": accession,
                    "id_method": method,
                    "via_table_s2_bridge": bridged,
                    "species": species,
                    "final_label": final_label,
                }
            )
            new_ind_lines.append(f"{sample_id:>30} {sex} {final_label}\n")

        write_atomic(Path(args.out), lambda h: h.writelines(new_ind_lines))

        def write_report(handle):
            writer = csv.DictWriter(
                handle,
                delimiter="\t",
                fieldnames=["sample_id", "accession", "id_method", "via_table_s2_bridge", "species", "final_label"],
            )
            writer.writeheader()
            writer.writerows(report_rows)

        write_atomic(Path(args.report), write_report)

    except (OSError, ValueError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(f"[labels] Table_S1.csv: {total_rows} data rows, {dup_count} duplicate Accession values (last wins)", file=sys.stderr)
    print(f"[labels] Table_S2.csv bridge: {len(table_s2)} SRA-id->Accession entries loaded, used for {bridge_used_count} samples", file=sys.stderr)
    print(
        f"[labels] total .ind samples: {len(ind_rows)} "
        f"(doubled-ID form: {method_counts.get('doubled', 0)}, plain-ID form: {method_counts.get('plain', 0)}), "
        f"matched Table_S1.csv: {n_matched}",
        file=sys.stderr,
    )
    print("[labels] final label (raw Group value) distribution:", file=sys.stderr)
    for label, count in label_counts.most_common():
        print(f"[labels]   {label}: {count}", file=sys.stderr)
    print("[labels] (species, group) combinations actually observed (for designing a standardization mapping later):", file=sys.stderr)
    for (species, group), count in species_group_counts.most_common():
        print(f"[labels]   species={species!r} group={group!r}: {count}", file=sys.stderr)
    print(f"[done] wrote {args.out}", file=sys.stderr)
    print(f"[done] wrote {args.report}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
