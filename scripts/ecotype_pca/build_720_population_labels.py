#!/usr/bin/env python3
"""
Attach real population labels (OrA-OrF wild-rice groups + cultivated
IND/AUS/ARO/TRJ/TEJ anchors) to asn720.6m.ind.

The panel's .ind file currently carries a single placeholder value,
"control", in its population column for all 720 samples (confirmed on the
server: `awk '{print $3}' asn720.6m.ind | sort -u` returns exactly one
value) -- this is docs/ECOTYPE_PCA_PANEL.md section 3.2 todo item 1.

TWO INDEPENDENT LABEL SOURCES, because asn720.6m.ind mixes two ID families
that are NOT resolvable the same way (confirmed against the real files,
not assumed):

  1. SRA-run-accession IDs ("ERR......" or "SRR......", 410 of 720
     samples) -> looked up directly by ID in asn720.pop.fam's IID column
     (2nd column). asn720.pop.fam's FID column (1st column) is used
     as-is as the label: OrA/OrB/OrC/OrD/OrE/OrF (wild-rice groups),
     OrADM, RAY, or the cultivated-rice anchor codes IND/AUS/ARO/TRJ/TEJ/
     ADM that this same panel also carries (asn720.pop.fam is a mixed
     wild+cultivated reference set, not wild-only, despite its name).
     asn720.pop.fam has 187 ERR-style + 533 SRR-style IDs, a superset of
     what's actually used in asn720.6m.ind.

  2. "IRIS_313-XXXXX_merged" / "B###_merged" IDs (310 of 720 samples,
     the "_merged" suffix marks accessions whose sequencing runs were
     merged into one BAM when this combined panel was built) -> these
     are the SAME physical 3K RGP accessions used in NB_final_snp.ind
     (see build_29m3k_population_labels.py). Stripping "_merged" recovers
     an ID directly resolvable in the same
     rice_line_metadata_20141029.xlsx metadata, via the exact same
     lookup/standardization logic -- reused here by importing
     build_29m3k_population_labels.py rather than duplicating it, so the
     two scripts can never silently drift into different label schemes
     for what is the same underlying data source. None of these 310
     samples exist in asn720.pop.fam under this ID form (confirmed: that
     file only contains ERR/SRR-style IIDs, never IRIS_313/B### form).

UNCONFIRMED LABEL MEANINGS -- flagged, not guessed at:
  - db/wild_rice_pangenome_README.txt (the file that was supposed to
    define OrA-OrF precisely, per docs/ECOTYPE_PCA_PANEL.md section 3.2
    todo item 1b) is EMPTY on the server. Whether OrA-OrF is the same
    numbering scheme as the "A pangenome reference of wild and cultivated
    rice" (Nature 2025) paper's Or-Ia/Or-Ib/Or-II/Or-IIIa/Or-IIIb groups
    is UNRESOLVED. Labels are passed through as OrA-OrF verbatim; do not
    assume a 1:1 mapping to the Nature 2025 paper's groups without an
    independent source.
  - "OrADM" is plausibly the wild-side equivalent of the cultivated
    "ADM" (admixed) label, by naming-pattern analogy alone -- not
    confirmed by any documentation found so far.
  - "RAY" (9 samples) has no evident meaning from any file found on the
    server. Passed through verbatim, unexplained.

SAFETY: same discipline as build_29m3k_population_labels.py -- hard
failure (non-zero exit, nothing written to --out) on duplicate IDs
(within asn720.6m.ind, and within asn720.pop.fam's IID column), a .ind
line that doesn't parse into exactly 3 fields, or a .pop.fam line that
doesn't parse into at least 2 fields. Atomic writes for the new .ind and
the coverage report.

Usage:
  python3 build_720_population_labels.py \\
    --ind /home/scratch/yinmt202607/db/6.7M_720/asn720.6m.ind \\
    --pop-fam /home/scratch/yinmt202607/db/asn720data/asn720.pop.fam \\
    --metadata-xlsx /home/scratch/yinmt202607/db/29M_3k/references/rice_line_metadata_20141029.xlsx \\
    --out /home/scratch/yinmt202607/db/6.7M_720/asn720.6m.labeled.ind \\
    --report /home/scratch/yinmt202607/db/6.7M_720/asn720.6m.label_report.tsv
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

# Same directory as this script when both are downloaded flat (matches how
# this repo's scripts are fetched onto the cluster, e.g. `gene/scripts/`).
import build_29m3k_population_labels as m3k

MERGED_SUFFIX = "_merged"
SRA_ID_RE = re.compile(r"^(ERR|SRR)\d+$")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--ind", required=True, help="path to asn720.6m.ind")
    parser.add_argument("--pop-fam", required=True, help="path to asn720.pop.fam")
    parser.add_argument("--metadata-xlsx", required=True, help="path to rice_line_metadata_20141029.xlsx")
    parser.add_argument("--out", required=True, help="path for the new .ind with real labels")
    parser.add_argument("--report", required=True, help="path for a per-sample coverage TSV")
    parser.add_argument(
        "--unmapped-label",
        default="UNK",
        help="population label for samples with no match in either source (default: UNK)",
    )
    return parser.parse_args(argv)


def validate_args(args: argparse.Namespace) -> None:
    for label, path in (("--ind", args.ind), ("--pop-fam", args.pop_fam), ("--metadata-xlsx", args.metadata_xlsx)):
        if not Path(path).is_file():
            raise FileNotFoundError(f"{label} file not found: {path}")


def load_pop_fam(path: Path) -> dict[str, str]:
    """Return IID -> FID (the population/group label), from a PLINK .fam file."""
    lookup: dict[str, str] = {}
    with path.open() as handle:
        for line_no, line in enumerate(handle, start=1):
            fields = line.split()
            if len(fields) < 2:
                raise ValueError(f"{path}:{line_no}: expected at least 2 fields (FID IID ...), got {line!r}")
            fid, iid = fields[0], fields[1]
            if iid in lookup:
                raise ValueError(f"duplicate IID in {path}: {iid!r}")
            lookup[iid] = fid
    return lookup


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


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        validate_args(args)
        pop_fam = load_pop_fam(Path(args.pop_fam))
        metadata = m3k.load_metadata(Path(args.metadata_xlsx))
        ind_rows = load_ind(Path(args.ind))

        report_rows = []
        new_ind_lines = []
        label_counts: Counter[str] = Counter()
        source_counts: Counter[str] = Counter()

        for sample_id, sex, _old_label in ind_rows:
            if sample_id.endswith(MERGED_SUFFIX):
                stripped = sample_id[: -len(MERGED_SUFFIX)]
                meta = metadata.get(stripped)
                if meta is None or meta[2] is None:
                    final_label = args.unmapped_label
                    source_sheet = meta[0] if meta is not None else ""
                    raw_group = meta[1] if meta is not None else ""
                    matched = False
                else:
                    source_sheet, raw_group, final_label = meta
                    matched = True
                source_type = "3k_metadata_via_merged_id"
            elif SRA_ID_RE.match(sample_id):
                fid = pop_fam.get(sample_id)
                if fid is None:
                    final_label = args.unmapped_label
                    matched = False
                else:
                    final_label = fid
                    matched = True
                source_sheet = ""
                raw_group = ""
                source_type = "asn720.pop.fam"
            else:
                final_label = args.unmapped_label
                matched = False
                source_sheet = ""
                raw_group = ""
                source_type = "unrecognized_id_format"

            source_counts[source_type] += 1
            if matched:
                source_counts[f"{source_type}_matched"] += 1
            label_counts[final_label] += 1
            report_rows.append(
                {
                    "sample_id": sample_id,
                    "source_type": source_type,
                    "matched": matched,
                    "source_sheet": source_sheet,
                    "raw_variety_group": raw_group if raw_group else "",
                    "final_label": final_label,
                }
            )
            new_ind_lines.append(f"{sample_id:>24} {sex} {final_label}\n")

        m3k.write_atomic(Path(args.out), lambda h: h.writelines(new_ind_lines))

        def write_report(handle):
            writer = csv.DictWriter(
                handle,
                delimiter="\t",
                fieldnames=["sample_id", "source_type", "matched", "source_sheet", "raw_variety_group", "final_label"],
            )
            writer.writeheader()
            writer.writerows(report_rows)

        m3k.write_atomic(Path(args.report), write_report)

    except (OSError, ValueError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(f"[labels] total samples: {len(ind_rows)}", file=sys.stderr)
    for source_type in ("asn720.pop.fam", "3k_metadata_via_merged_id", "unrecognized_id_format"):
        total = source_counts.get(source_type, 0)
        if total == 0:
            continue
        matched = source_counts.get(f"{source_type}_matched", 0)
        print(f"[labels]   {source_type}: {matched}/{total} matched", file=sys.stderr)
    print("[labels] final label distribution:", file=sys.stderr)
    for label, count in label_counts.most_common():
        print(f"[labels]   {label}: {count}", file=sys.stderr)
    print(f"[done] wrote {args.out}", file=sys.stderr)
    print(f"[done] wrote {args.report}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
