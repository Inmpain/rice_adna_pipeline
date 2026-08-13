#!/usr/bin/env python3
"""
Attach real population labels to civan_snp.ind from Civáň et al. 2019's
own Table_S1.csv sample metadata.

The panel's .ind file currently carries a single placeholder value,
"???", in its population column for all 1,056 samples -- this is
docs/ECOTYPE_PCA_PANEL.md section 3.2 todo item 2c.

ID FORMAT: every civan_snp.ind sample ID is the underlying accession
name duplicated around a literal underscore -- "B006" becomes
"B006_B006", "IRIS_313-9986" (itself containing an underscore) becomes
"IRIS_313-9986_IRIS_313-9986". This is plink2's default FID_IID
behavior when converting a VCF whose samples have no separate
family/individual ID (see docs/ECOTYPE_PCA_PANEL.md 1.3/5.6). Recovering
the real accession therefore CANNOT be done by splitting on the first
or last underscore (that breaks for accessions like "IRIS_313-9986"
which contain underscores themselves) -- it must find the exact
midpoint: the string has odd length 2n+1, character n is '_', and the
two n-character halves are identical. recover_accession() below does
exactly that and returns None (hard error) if a sample ID doesn't fit
this pattern, since every sample so far has been confirmed to follow it.

WHY THIS SCRIPT DOES NOT STANDARDIZE THE LABEL, UNLIKE THE OTHER TWO
build_*_population_labels.py SCRIPTS: Table_S1.csv's header is messy
(quoted fields containing literal newlines break naive line-based
reading -- must use Python's csv module) and this script's author could
not get a clean, exhaustive list of the "Group" column's actual values
from the server before writing this (unlike NB_final_snp.ind's 8-value
xlsx column, which was fully scanned first). The raw "Group" value is
therefore used AS THE LABEL VERBATIM, and "Species" is carried in the
report alongside it for context. Once a real run's [labels] stderr
summary shows the true value distribution, a standardization mapping
(if one is even needed -- Civáň's Group column may already use clean
values like "indica"/"japonica"/"aus"/"aromatic") can be added as a
followup, the same way NB_final_snp.ind's mapping was built from an
observed value list rather than guessed in advance.

Table_S1.csv is also known (docs/ECOTYPE_PCA_PANEL.md 1.3) to have 1,063
data rows against a paper-stated 1,056 samples -- 7 more than expected,
cause unconfirmed. Rather than hard-failing on a duplicate Accession (as
the other two scripts do, where duplicates would be a genuine anomaly),
this script tolerates duplicate Accession values in Table_S1.csv
(last-row-wins) and reports the duplicate count, since this file is
already a known, unexplained edge case -- a hard failure here would
block ever seeing the real numbers needed to explain it.

Usage:
  python3 build_civan_population_labels.py \\
    --ind /home/scratch/yinmt202607/db/paper1/civan_snp.ind \\
    --table-s1 /home/scratch/yinmt202607/db/paper1/Table_S1.csv \\
    --out /home/scratch/yinmt202607/db/paper1/civan_snp.labeled.ind \\
    --report /home/scratch/yinmt202607/db/paper1/civan_snp.label_report.tsv
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
import tempfile
from collections import Counter
from pathlib import Path


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--ind", required=True, help="path to civan_snp.ind")
    parser.add_argument("--table-s1", required=True, help="path to Civáň et al. 2019 Table_S1.csv")
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


def recover_accession(sample_id: str) -> str | None:
    n = len(sample_id)
    if n % 2 == 0:
        return None
    mid = n // 2
    if sample_id[mid] != "_":
        return None
    first_half, second_half = sample_id[:mid], sample_id[mid + 1 :]
    if first_half != second_half:
        return None
    return first_half


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
        ind_rows = load_ind(Path(args.ind))

        report_rows = []
        new_ind_lines = []
        label_counts: Counter[str] = Counter()
        species_group_counts: Counter[tuple[str, str]] = Counter()
        n_recovered = 0
        n_matched = 0

        for sample_id, sex, _old_label in ind_rows:
            accession = recover_accession(sample_id)
            if accession is None:
                raise ValueError(
                    f"sample ID {sample_id!r} does not fit the expected "
                    "'ACCESSION_ACCESSION' self-duplicated pattern -- every "
                    "sample seen so far has fit it, so this is treated as a "
                    "hard error rather than silently guessed at."
                )
            n_recovered += 1
            meta = table_s1.get(accession)
            if meta is None:
                final_label = args.unmapped_label
                species, group = "", ""
            else:
                species, group = meta
                final_label = group if group else args.unmapped_label
                n_matched += 1
                species_group_counts[(species, group)] += 1
            label_counts[final_label] += 1
            report_rows.append(
                {
                    "sample_id": sample_id,
                    "accession": accession,
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
                fieldnames=["sample_id", "accession", "species", "final_label"],
            )
            writer.writeheader()
            writer.writerows(report_rows)

        write_atomic(Path(args.report), write_report)

    except (OSError, ValueError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(f"[labels] Table_S1.csv: {total_rows} data rows, {dup_count} duplicate Accession values (last wins)", file=sys.stderr)
    print(f"[labels] total .ind samples: {len(ind_rows)}, accession recovered: {n_recovered}, matched Table_S1.csv: {n_matched}", file=sys.stderr)
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
