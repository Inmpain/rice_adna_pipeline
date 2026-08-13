#!/usr/bin/env python3
"""Split genus-wide best-hit reads into a panel-target set and other Oryza.

The existing best-hit run classified every retained read to its best Oryza
species and wrote that assignment to ``oryza_filter.decisions.tsv.gz``.  This
script reuses those decisions; it does not rerun competitive mapping.

The default target is deliberately the Oryza rufipogon species complex used by
the current population panels: O. rufipogon (4529), O. sativa (4530), and
O. nivara (4536).  It is not labelled "all AA-genome Oryza", because the
current panels do not represent all AA-genome species.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import sys
import tempfile
from collections import Counter
from contextlib import ExitStack
from pathlib import Path
from typing import TextIO


DEFAULT_TARGET_TAXIDS = (4529, 4530, 4536)
REQUIRED_DECISION_COLUMNS = {
    "read_name",
    "best_oryza_taxid",
    "best_oryza_name",
    "decision",
}


def parse_taxids(value: str) -> frozenset[int]:
    fields = value.replace(",", " ").split()
    if not fields:
        raise argparse.ArgumentTypeError("target taxid list must not be empty")
    try:
        taxids = frozenset(int(field) for field in fields)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "target taxids must be comma- or space-separated integers"
        ) from exc
    if any(taxid <= 0 for taxid in taxids):
        raise argparse.ArgumentTypeError("target taxids must be positive integers")
    return taxids


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument(
        "--besthit-fastq",
        required=True,
        help="genus-wide <sample>.besthit_oryza.fastq.gz",
    )
    parser.add_argument(
        "--decisions",
        required=True,
        help="matching <sample>.oryza_filter.decisions.tsv.gz",
    )
    parser.add_argument("--outdir", required=True)
    parser.add_argument(
        "--target-taxids",
        type=parse_taxids,
        default=frozenset(DEFAULT_TARGET_TAXIDS),
        help="target species taxids (default: 4529,4530,4536)",
    )
    parser.add_argument(
        "--target-label",
        default="target_orsc",
        help="filename label for target reads (default: target_orsc)",
    )
    parser.add_argument(
        "--other-label",
        default="other_oryza",
        help="filename label for non-target retained reads",
    )
    return parser.parse_args(argv)


def open_text(path: Path, mode: str) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, mode + "t", newline="")
    return path.open(mode, newline="")


def load_kept_assignments(
    decisions_path: Path,
) -> tuple[dict[str, tuple[int, str]], Counter[tuple[int, str]]]:
    assignments: dict[str, tuple[int, str]] = {}
    species_counts: Counter[tuple[int, str]] = Counter()
    with open_text(decisions_path, "r") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = REQUIRED_DECISION_COLUMNS - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"decision table is missing columns: {', '.join(sorted(missing))}"
            )
        for line_number, row in enumerate(reader, start=2):
            if row["decision"] != "KEEP":
                continue
            read_name = row["read_name"]
            if read_name in assignments:
                raise ValueError(
                    f"duplicate KEEP read_name at decision line {line_number}: {read_name}"
                )
            raw_taxid = row["best_oryza_taxid"]
            if raw_taxid in {"", "NA"}:
                raise ValueError(
                    f"KEEP row has no best_oryza_taxid at line {line_number}: {read_name}"
                )
            try:
                taxid = int(raw_taxid)
            except ValueError as exc:
                raise ValueError(
                    f"invalid best_oryza_taxid at line {line_number}: {raw_taxid}"
                ) from exc
            species_name = row["best_oryza_name"]
            assignments[read_name] = (taxid, species_name)
            species_counts[(taxid, species_name)] += 1
    return assignments, species_counts


def fastq_name(header: str, record_number: int) -> str:
    if not header.startswith("@"):
        raise ValueError(
            f"FASTQ record {record_number} header does not start with '@': {header.rstrip()}"
        )
    name = header[1:].split(maxsplit=1)[0]
    if not name:
        raise ValueError(f"FASTQ record {record_number} has an empty read name")
    return name


def write_outputs(args: argparse.Namespace) -> tuple[Path, Path, Path, Path]:
    besthit_path = Path(args.besthit_fastq)
    decisions_path = Path(args.decisions)
    if not besthit_path.is_file():
        raise FileNotFoundError(f"best-hit FASTQ not found: {besthit_path}")
    if not decisions_path.is_file():
        raise FileNotFoundError(f"decision table not found: {decisions_path}")

    assignments, decision_species_counts = load_kept_assignments(decisions_path)
    if not assignments:
        raise ValueError("decision table contains no KEEP reads")

    outdir = Path(args.outdir)
    target_path = outdir / f"{args.sample}.{args.target_label}.fastq.gz"
    other_path = outdir / f"{args.sample}.{args.other_label}.fastq.gz"
    summary_path = outdir / f"{args.sample}.taxonomic_tiers.summary.tsv"
    species_path = outdir / f"{args.sample}.taxonomic_tiers.by_species.tsv"

    outdir.mkdir(parents=True, exist_ok=True)
    temp_paths: dict[Path, Path] = {}
    try:
        with ExitStack() as stack:
            handles: dict[Path, TextIO] = {}
            for final_path in (target_path, other_path, summary_path, species_path):
                suffix = ".tmp.gz" if final_path.suffix == ".gz" else ".tmp"
                fd, tmp_name = tempfile.mkstemp(
                    prefix=f".{final_path.name}.", suffix=suffix, dir=outdir
                )
                os.close(fd)
                tmp_path = Path(tmp_name)
                temp_paths[final_path] = tmp_path
                handles[final_path] = stack.enter_context(open_text(tmp_path, "w"))

            observed: set[str] = set()
            category_counts: Counter[str] = Counter()
            observed_species_counts: Counter[tuple[str, int, str]] = Counter()
            with open_text(besthit_path, "r") as fastq:
                record_number = 0
                while True:
                    header = fastq.readline()
                    if not header:
                        break
                    sequence = fastq.readline()
                    plus = fastq.readline()
                    quality = fastq.readline()
                    record_number += 1
                    if not sequence or not plus or not quality:
                        raise ValueError(f"truncated FASTQ at record {record_number}")
                    if not plus.startswith("+"):
                        raise ValueError(
                            f"FASTQ record {record_number} third line does not start with '+'"
                        )
                    if len(sequence.rstrip("\r\n")) != len(quality.rstrip("\r\n")):
                        raise ValueError(
                            f"sequence/quality length mismatch at FASTQ record {record_number}"
                        )
                    read_name = fastq_name(header, record_number)
                    if read_name in observed:
                        raise ValueError(f"duplicate read name in best-hit FASTQ: {read_name}")
                    observed.add(read_name)
                    if read_name not in assignments:
                        raise ValueError(
                            f"best-hit FASTQ read has no KEEP assignment: {read_name}"
                        )
                    taxid, species_name = assignments[read_name]
                    category = "target" if taxid in args.target_taxids else "other_oryza"
                    output = handles[target_path] if category == "target" else handles[other_path]
                    output.write(header)
                    output.write(sequence)
                    output.write(plus)
                    output.write(quality)
                    category_counts[category] += 1
                    observed_species_counts[(category, taxid, species_name)] += 1

            missing_fastq = set(assignments) - observed
            if missing_fastq:
                example = min(missing_fastq)
                raise ValueError(
                    f"{len(missing_fastq)} KEEP decisions are absent from the best-hit FASTQ; "
                    f"example: {example}"
                )
            if sum(decision_species_counts.values()) != len(observed):
                raise AssertionError("internal KEEP/species count inconsistency")

            total = len(observed)
            summary = handles[summary_path]
            summary.write("sample\ttarget_label\ttarget_taxids\tbesthit_kept_reads\t")
            summary.write("target_reads\tother_oryza_reads\ttarget_pct_of_kept\n")
            target_count = category_counts["target"]
            other_count = category_counts["other_oryza"]
            target_pct = 100.0 * target_count / total if total else 0.0
            summary.write(
                f"{args.sample}\t{args.target_label}\t"
                f"{','.join(str(x) for x in sorted(args.target_taxids))}\t{total}\t"
                f"{target_count}\t{other_count}\t{target_pct:.2f}\n"
            )

            by_species = handles[species_path]
            by_species.write("sample\tcategory\ttaxid\tspecies_name\treads\tpct_of_kept\n")
            for (category, taxid, species_name), count in sorted(
                observed_species_counts.items(), key=lambda item: (-item[1], item[0])
            ):
                pct = 100.0 * count / total if total else 0.0
                by_species.write(
                    f"{args.sample}\t{category}\t{taxid}\t{species_name}\t{count}\t{pct:.2f}\n"
                )

        for final_path, tmp_path in temp_paths.items():
            os.replace(tmp_path, final_path)
    except BaseException:
        for tmp_path in temp_paths.values():
            tmp_path.unlink(missing_ok=True)
        raise
    return target_path, other_path, summary_path, species_path


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        outputs = write_outputs(args)
    except (OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    print(f"[done] sample={args.sample} target_taxids={','.join(map(str, sorted(args.target_taxids)))}")
    for output in outputs:
        print(f"[done] wrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
