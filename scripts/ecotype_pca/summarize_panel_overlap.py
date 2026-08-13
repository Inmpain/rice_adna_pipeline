#!/usr/bin/env python3
"""Count panel SNPs observed in one sparse ancient-rice BAM.

The BAM is scanned once. Duplicate, unmapped, secondary, supplementary and
QC-failed alignments are excluded. Base quality is applied while two mapping-
quality tracks (default Q0 and Q20) are retained side by side. Each EIGENSTRAT
``.snp`` panel is then streamed without loading its millions of rows in memory.

"covered" means at least one eligible base overlaps the panel coordinate.
"allele_supported" means at least one eligible base equals either panel allele.
The distinction prevents coverage from being mistaken for a callable panel SNP.
"callable_tv" applies the same transversion-only restriction as the main
pseudo-haploid analysis; "callable_all" includes transitions as well.

A ``[qc]`` line is always printed to stderr breaking down every pileup entry
seen in the BAM by why it was excluded (read flags, missing query position,
sub-threshold base quality, sub-threshold mapping quality) versus how many
bases actually passed each mapping-quality track. Pass ``--qc-out`` to also
write this breakdown as a one-row TSV, so filter-driven data loss on a thin
ancient-DNA sample is visible instead of silently dropped.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

try:
    import pysam
except ImportError:  # pragma: no cover - exercised on cluster, not stdlib CI
    pysam = None


TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
BASES = frozenset("ACGT")


@dataclass(frozen=True)
class Panel:
    label: str
    path: Path


def parse_panel(value: str) -> Panel:
    if "=" not in value:
        raise argparse.ArgumentTypeError("panel must be LABEL=/path/to/panel.snp")
    label, raw_path = value.split("=", 1)
    if not label or not re.fullmatch(r"[A-Za-z0-9_.-]+", label):
        raise argparse.ArgumentTypeError(
            "panel label must contain only letters, digits, dot, underscore, or dash"
        )
    if not raw_path:
        raise argparse.ArgumentTypeError("panel path must not be empty")
    return Panel(label, Path(raw_path))


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--bam", required=True)
    parser.add_argument(
        "--panel",
        action="append",
        type=parse_panel,
        required=True,
        help="repeatable LABEL=/path/to/panel.snp",
    )
    parser.add_argument("--out", required=True)
    parser.add_argument(
        "--qc-out",
        default=None,
        help="optional path to write per-sample BAM-level quality-filter counts as TSV",
    )
    parser.add_argument("--min-baseq", type=int, default=20)
    parser.add_argument("--low-mapq", type=int, default=0)
    parser.add_argument("--high-mapq", type=int, default=20)
    return parser.parse_args(argv)


def validate_args(args: argparse.Namespace) -> None:
    if pysam is None:
        raise RuntimeError(
            "pysam is required (activate the same Python environment used for "
            "pseudo_haploid_call.py)"
        )
    bam = Path(args.bam)
    if not bam.is_file():
        raise FileNotFoundError(f"BAM not found: {bam}")
    if args.min_baseq < 0 or args.low_mapq < 0 or args.high_mapq < 0:
        raise ValueError("quality thresholds must be non-negative")
    if args.low_mapq > args.high_mapq:
        raise ValueError("--low-mapq must be <= --high-mapq")
    labels = [panel.label for panel in args.panel]
    if len(labels) != len(set(labels)):
        raise ValueError("panel labels must be unique")
    for panel in args.panel:
        if not panel.path.is_file():
            raise FileNotFoundError(f"panel SNP file not found: {panel.path}")


def normalize_contig(value: str, bam_references: frozenset[str]) -> str | None:
    candidates = [value]
    stripped = value.removeprefix("chr").removeprefix("Chr")
    if stripped.isdigit():
        number = int(stripped)
        candidates.extend(
            [str(number), f"{number:02d}", f"chr{number}", f"chr{number:02d}", f"Chr{number}"]
        )
    for candidate in candidates:
        if candidate in bam_references:
            return candidate
    return None


def build_coverage_index(
    bam_path: Path,
    min_baseq: int,
    low_mapq: int,
    high_mapq: int,
) -> tuple[dict[tuple[str, int], tuple[frozenset[str], frozenset[str]]], Counter[str]]:
    coverage: dict[tuple[str, int], tuple[frozenset[str], frozenset[str]]] = {}
    qc: Counter[str] = Counter()
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        qc["bam_references"] = len(bam.references)
        for contig in bam.references:
            for column in bam.pileup(
                contig,
                stepper="all",
                min_base_quality=0,
                min_mapping_quality=0,
                ignore_overlaps=False,
                ignore_orphans=False,
                truncate=False,
                max_depth=1_000_000,
                # pysam's own default flag_filter (0x704: unmapped, secondary,
                # qcfail, duplicate) silently drops those reads before this
                # loop ever sees them, making the excluded_flag_* QC counter
                # below always read 0 regardless of real duplicate content.
                # Disable it so the explicit flag checks in this loop are the
                # only filtering step, and the QC counts reflect reality.
                flag_filter=0,
            ):
                low_bases: set[str] = set()
                high_bases: set[str] = set()
                for pileup_read in column.pileups:
                    qc["pileup_entries_total"] += 1
                    if pileup_read.is_del or pileup_read.is_refskip:
                        qc["excluded_del_or_refskip"] += 1
                        continue
                    aln = pileup_read.alignment
                    if (
                        aln.is_unmapped
                        or aln.is_duplicate
                        or aln.is_secondary
                        or aln.is_supplementary
                        or aln.is_qcfail
                    ):
                        qc["excluded_flag_dup_secondary_supp_qcfail_unmapped"] += 1
                        continue
                    query_position = pileup_read.query_position
                    if query_position is None or aln.query_sequence is None:
                        qc["excluded_no_query_position"] += 1
                        continue
                    qualities = aln.query_qualities
                    if qualities is None or qualities[query_position] < min_baseq:
                        qc["excluded_low_baseq"] += 1
                        continue
                    base = aln.query_sequence[query_position].upper()
                    if base not in BASES:
                        qc["excluded_non_acgt_base"] += 1
                        continue
                    if aln.mapping_quality < low_mapq:
                        qc["excluded_low_mapq_below_low"] += 1
                        continue
                    low_bases.add(base)
                    qc["bases_pass_low_mapq"] += 1
                    if aln.mapping_quality >= high_mapq:
                        high_bases.add(base)
                        qc["bases_pass_high_mapq"] += 1
                    else:
                        qc["bases_mid_mapq_low_only"] += 1
                if low_bases:
                    coverage[(contig, column.reference_pos + 1)] = (
                        frozenset(low_bases),
                        frozenset(high_bases),
                    )
        qc["covered_genome_positions_low"] = len(coverage)
        qc["covered_genome_positions_high"] = sum(
            bool(high) for _, high in coverage.values()
        )
    return coverage, qc


def empty_panel_counts() -> Counter[str]:
    return Counter(
        {
            "total_panel_snps": 0,
            "valid_biallelic_snps": 0,
            "transversion_snps": 0,
            "invalid_panel_rows": 0,
            "unresolved_contig_rows": 0,
            "covered_low": 0,
            "covered_high": 0,
            "allele_supported_low": 0,
            "allele_supported_high": 0,
            "allele_mismatch_only_low": 0,
            "allele_mismatch_only_high": 0,
            "callable_all_low": 0,
            "callable_all_high": 0,
            "callable_tv_low": 0,
            "callable_tv_high": 0,
        }
    )


def count_panel(
    panel: Panel,
    coverage: dict[tuple[str, int], tuple[frozenset[str], frozenset[str]]],
    bam_references: frozenset[str],
) -> tuple[Counter[str], set[str], set[str]]:
    counts = empty_panel_counts()
    chromosomes_low: set[str] = set()
    chromosomes_high: set[str] = set()
    with panel.path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            counts["total_panel_snps"] += 1
            fields = line.split()
            if len(fields) < 6:
                counts["invalid_panel_rows"] += 1
                continue
            try:
                position = int(fields[3])
            except ValueError:
                counts["invalid_panel_rows"] += 1
                continue
            allele1, allele2 = fields[4].upper(), fields[5].upper()
            if allele1 not in BASES or allele2 not in BASES or allele1 == allele2:
                counts["invalid_panel_rows"] += 1
                continue
            counts["valid_biallelic_snps"] += 1
            is_tv = (allele1, allele2) not in TRANSITIONS
            if is_tv:
                counts["transversion_snps"] += 1
            contig = normalize_contig(fields[1], bam_references)
            if contig is None:
                counts["unresolved_contig_rows"] += 1
                continue
            bases = coverage.get((contig, position))
            if bases is None:
                continue
            low_bases, high_bases = bases
            alleles = {allele1, allele2}
            counts["covered_low"] += 1
            if low_bases & alleles:
                counts["allele_supported_low"] += 1
                counts["callable_all_low"] += 1
                chromosomes_low.add(contig)
                if is_tv:
                    counts["callable_tv_low"] += 1
            else:
                counts["allele_mismatch_only_low"] += 1
            if high_bases:
                counts["covered_high"] += 1
                if high_bases & alleles:
                    counts["allele_supported_high"] += 1
                    counts["callable_all_high"] += 1
                    chromosomes_high.add(contig)
                    if is_tv:
                        counts["callable_tv_high"] += 1
                else:
                    counts["allele_mismatch_only_high"] += 1
    return counts, chromosomes_low, chromosomes_high


def write_summary(
    args: argparse.Namespace,
    coverage: dict[tuple[str, int], tuple[frozenset[str], frozenset[str]]],
    qc: Counter[str],
) -> None:
    # Recover all BAM reference names directly from the header; chromosomes
    # with no retained bases are absent from the sparse coverage dictionary.
    with pysam.AlignmentFile(args.bam, "rb") as bam:
        bam_references = frozenset(bam.references)

    rows = []
    for panel in args.panel:
        counts, chrom_low, chrom_high = count_panel(panel, coverage, bam_references)
        rows.append((panel, counts, chrom_low, chrom_high))

    output = Path(args.out)
    output.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=f".{output.name}.", suffix=".tmp", dir=output.parent)
    os.close(fd)
    tmp_path = Path(tmp_name)
    fields = [
        "sample", "panel", "panel_snp_file", "min_baseq", "low_mapq", "high_mapq",
        "genome_positions_low", "genome_positions_high", "total_panel_snps",
        "valid_biallelic_snps", "transversion_snps", "invalid_panel_rows",
        "unresolved_contig_rows", "covered_low", "covered_high",
        "allele_supported_low", "allele_supported_high", "allele_mismatch_only_low",
        "allele_mismatch_only_high", "callable_all_low", "callable_all_high",
        "callable_tv_low", "callable_tv_high",
        "chromosomes_callable_low", "chromosomes_callable_high",
    ]
    try:
        with tmp_path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
            writer.writeheader()
            for panel, counts, chrom_low, chrom_high in rows:
                row = {
                    "sample": args.sample,
                    "panel": panel.label,
                    "panel_snp_file": str(panel.path),
                    "min_baseq": args.min_baseq,
                    "low_mapq": args.low_mapq,
                    "high_mapq": args.high_mapq,
                    "genome_positions_low": qc["covered_genome_positions_low"],
                    "genome_positions_high": qc["covered_genome_positions_high"],
                    **{field: counts[field] for field in counts},
                    "chromosomes_callable_low": len(chrom_low),
                    "chromosomes_callable_high": len(chrom_high),
                }
                writer.writerow(row)
        os.replace(tmp_path, output)
    except BaseException:
        tmp_path.unlink(missing_ok=True)
        raise


QC_FIELDS = [
    "sample", "min_baseq", "low_mapq", "high_mapq", "bam_references",
    "pileup_entries_total", "excluded_del_or_refskip",
    "excluded_flag_dup_secondary_supp_qcfail_unmapped",
    "excluded_no_query_position", "excluded_low_baseq", "excluded_non_acgt_base",
    "excluded_low_mapq_below_low", "bases_mid_mapq_low_only", "bases_pass_low_mapq",
    "bases_pass_high_mapq", "covered_genome_positions_low", "covered_genome_positions_high",
]


def qc_summary_line(sample: str, qc: Counter[str]) -> str:
    return (
        f"[qc] {sample}: pileup_entries={qc['pileup_entries_total']} "
        f"excluded_flags={qc['excluded_flag_dup_secondary_supp_qcfail_unmapped']} "
        f"excluded_low_baseq={qc['excluded_low_baseq']} "
        f"excluded_low_mapq={qc['excluded_low_mapq_below_low']} "
        f"bases_pass_low_mapq={qc['bases_pass_low_mapq']} "
        f"bases_pass_high_mapq={qc['bases_pass_high_mapq']}"
    )


def write_qc_summary(args: argparse.Namespace, qc: Counter[str]) -> None:
    output = Path(args.qc_out)
    output.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=f".{output.name}.", suffix=".tmp", dir=output.parent)
    os.close(fd)
    tmp_path = Path(tmp_name)
    try:
        with tmp_path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, delimiter="\t", fieldnames=QC_FIELDS)
            writer.writeheader()
            row = {
                "sample": args.sample,
                "min_baseq": args.min_baseq,
                "low_mapq": args.low_mapq,
                "high_mapq": args.high_mapq,
            }
            row.update({field: qc[field] for field in QC_FIELDS if field not in row})
            writer.writerow(row)
        os.replace(tmp_path, output)
    except BaseException:
        tmp_path.unlink(missing_ok=True)
        raise


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        validate_args(args)
        print(f"[overlap] indexing {args.bam}", file=sys.stderr)
        coverage, qc = build_coverage_index(
            Path(args.bam), args.min_baseq, args.low_mapq, args.high_mapq
        )
        write_summary(args, coverage, qc)
        print(qc_summary_line(args.sample, qc), file=sys.stderr)
        if args.qc_out:
            write_qc_summary(args, qc)
            print(f"[qc] wrote {args.qc_out}", file=sys.stderr)
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    print(f"[done] wrote {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
