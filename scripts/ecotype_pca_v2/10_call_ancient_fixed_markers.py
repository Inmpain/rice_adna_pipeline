#!/usr/bin/env python3
"""Call one deterministic pseudo-haploid allele at each frozen marker."""

from __future__ import annotations

import argparse
import csv
import random
import sys
from collections import defaultdict
from pathlib import Path

from fixed_projection_lib import POOLED_LIBRARY_TYPE, iter_snp, refuse_existing, stable_int_seed, write_tsv
from lib_ecotype_v2 import load_config


TRANSITIONS = {frozenset(("A", "G")), frozenset(("C", "T"))}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--fixed-snp", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--panel", choices=("A", "B", "C"), required=True)
    parser.add_argument("--library-type", choices=(POOLED_LIBRARY_TYPE,), required=True)
    parser.add_argument("--track", choices=("TV", "ALL"), required=True)
    parser.add_argument("--contig-format", default="chr%02d")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def site_seed(sample: str, panel: str, contig: str, pos: int) -> int:
    """Track/library intentionally absent so TV and ALL share draws."""
    return stable_int_seed(sample, panel, contig, pos)


def main() -> int:
    args = parse_args()
    try:
        import pysam
    except ImportError:
        print("FATAL: pysam is required", file=sys.stderr)
        return 2
    cfg = load_config(args.config)
    mapq, baseq = int(cfg["ancient"]["mapq"]), int(cfg["ancient"]["baseq"])
    bam_path = Path(args.bam)
    if not bam_path.is_file():
        print(f"FATAL: BAM not found: {bam_path}", file=sys.stderr)
        return 3
    records = list(iter_snp(args.fixed_snp))
    if not records or any(record["ref"] is None for record in records):
        print("FATAL: fixed .snp must be nonempty and have six columns", file=sys.stderr)
        return 3
    if args.track == "TV":
        transitions = [record["id"] for record in records if frozenset((record["ref"], record["alt"])) in TRANSITIONS]
        if transitions:
            print(f"FATAL: TV fixed list contains transition SNPs: {transitions[:10]}", file=sys.stderr)
            return 3

    output_dir = Path(args.out_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = output_dir / f"{args.sample}.{args.panel}.{args.library_type}.{args.track}"
    calls_path = Path(str(prefix) + ".calls.txt")
    sites_path = Path(str(prefix) + ".call_sites.tsv")
    report_path = Path(str(prefix) + ".call_report.tsv")
    try:
        refuse_existing([calls_path, sites_path, report_path], args.overwrite)
        by_contig = defaultdict(set)
        contig_for = {}
        for record in records:
            try:
                contig = args.contig_format % int(record["chrom"].lstrip("0") or "0")
            except (TypeError, ValueError) as exc:
                raise ValueError(f"unsupported chromosome {record['chrom']!r}") from exc
            by_contig[contig].add(record["pos"])
            contig_for[(record["chrom"], record["pos"])] = contig

        chosen = {}
        drawn_count = defaultdict(int)
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            references = set(bam.references)
            absent = sorted(set(by_contig) - references)
            if len(absent) == len(by_contig):
                raise ValueError(f"no formatted panel contigs occur in BAM; tried {absent[:5]}")
            if bam.mapped == 0:
                raise ValueError("BAM has zero mapped reads")
            for contig, targets in by_contig.items():
                if contig not in references:
                    continue
                for column in bam.pileup(contig, min_mapping_quality=mapq, stepper="samtools", ignore_overlaps=False, ignore_orphans=False):
                    pos = column.reference_pos + 1
                    if pos not in targets:
                        continue
                    rng = random.Random(site_seed(args.sample, args.panel, contig, pos))
                    selected = None
                    n = 0
                    for pileup_read in column.pileups:
                        if pileup_read.is_del or pileup_read.is_refskip:
                            continue
                        alignment = pileup_read.alignment
                        if alignment.is_duplicate or alignment.is_secondary or alignment.is_supplementary or alignment.is_qcfail:
                            continue
                        if alignment.mapping_quality < mapq:
                            continue
                        query_pos = pileup_read.query_position
                        if query_pos is None or alignment.query_qualities is None:
                            continue
                        if alignment.query_qualities[query_pos] < baseq:
                            continue
                        base = alignment.query_sequence[query_pos].upper()
                        n += 1
                        if rng.randrange(n) == 0:
                            selected = base
                    if selected is not None:
                        chosen[(contig, pos)] = selected
                        drawn_count[(contig, pos)] = n

        site_rows, calls = [], []
        counts = defaultdict(int)
        for record in records:
            contig = contig_for[(record["chrom"], record["pos"])]
            base = chosen.get((contig, record["pos"]))
            if base is None:
                call, status = "9", "NO_COVERAGE"
            elif base == record["ref"]:
                call, status = "2", "CALLED_REF"
            elif base == record["alt"]:
                call, status = "0", "CALLED_ALT"
            else:
                call, status = "9", "ALLELE_MISMATCH"
            calls.append(call)
            counts[status] += 1
            site_rows.append({
                "snp_id": record["id"], "chrom": record["chrom"], "position": record["pos"],
                "ref": record["ref"], "alt": record["alt"], "call": call, "status": status,
                "qualifying_read_n": drawn_count.get((contig, record["pos"]), 0),
            })
        calls_path.write_text("\n".join(calls) + "\n")
        write_tsv(sites_path, site_rows, list(site_rows[0]))
        callable_n = counts["CALLED_REF"] + counts["CALLED_ALT"]
        drawn_n = callable_n + counts["ALLELE_MISMATCH"]
        metrics = [
            ("sample", args.sample), ("panel", args.panel), ("library_type", args.library_type),
            ("track", args.track), ("fixed_marker_n", len(records)), ("mapq", mapq), ("baseq", baseq),
            ("called", callable_n), ("no_coverage", counts["NO_COVERAGE"]),
            ("allele_mismatch", counts["ALLELE_MISMATCH"]),
            ("callable_fraction", f"{callable_n / len(records):.10f}"),
            ("allele_match_rate_among_covered", f"{(callable_n / drawn_n) if drawn_n else 0.0:.10f}"),
            ("seed_contract", "sha256(sample:panel:contig:position);track_excluded"),
            ("bam_interpretation", "pooled_mixed_capture_plus_shotgun"),
        ]
        with open(report_path, "w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(("metric", "value"))
            writer.writerows(metrics)
    except (OSError, ValueError) as exc:
        print(f"FATAL: {exc}", file=sys.stderr)
        return 3
    print(f"PASS: {args.sample} {args.track}: {callable_n}/{len(records)} fixed markers callable")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
