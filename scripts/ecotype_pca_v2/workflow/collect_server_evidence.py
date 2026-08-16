#!/usr/bin/env python3
"""Collect small, read-only evidence needed to resolve current v2 blockers.

This deliberately records facts instead of making data decisions: exact raw
versus filtered sample IDs/labels and paired-read flag counts for each ancient
BAM.  It never reads or copies genotype matrices into its output bundle.
"""

import argparse
import csv
import json
import re
import subprocess
from pathlib import Path


def read_ind(path):
    rows = {}
    with open(path) as handle:
        for lineno, line in enumerate(handle, 1):
            fields = line.split()
            if len(fields) < 3:
                raise ValueError(f"{path}:{lineno}: expected ID sex label")
            if fields[0] in rows:
                raise ValueError(f"{path}:{lineno}: duplicate sample ID {fields[0]}")
            rows[fields[0]] = {"sex": fields[1], "label": fields[2]}
    return rows


def parse_flagstat(text):
    """Return total/paired/proper-pair counts from one samtools flagstat pass."""
    counts = {}
    patterns = {
        "records": "in total",
        "paired_flag": "paired in sequencing",
        "proper_pair_flag": "properly paired",
    }
    for line in text.splitlines():
        match = re.match(r"^(\d+) \+ (\d+) (.+)$", line.strip())
        if not match:
            continue
        combined = int(match.group(1)) + int(match.group(2))
        label = match.group(3)
        for key, prefix in patterns.items():
            if label.startswith(prefix):
                counts[key] = combined
    missing = sorted(set(patterns) - set(counts))
    if missing:
        raise ValueError(f"samtools flagstat output missing fields: {missing}")
    return counts


def samtools_flagstat(bam):
    proc = subprocess.run(
        ["samtools", "flagstat", str(bam)],
        capture_output=True, text=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"samtools failed for {bam}: {proc.stderr.strip()}")
    return parse_flagstat(proc.stdout)


def main():
    import yaml

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()

    cfg = yaml.safe_load(open(args.config))
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=False)
    evidence = {"panels": {}, "ancient_bams": []}

    diff_rows = []
    for panel_key in ("panel_A_3k", "panel_B_720", "panel_C_civan"):
        info = cfg["inputs"][panel_key]
        root = Path(info["dir"])
        prefix = info["prefix"]
        suffix = info["filtered_suffix"]
        raw_path = root / f"{prefix}.ind"
        filtered_path = root / f"{prefix}{suffix}.ind"
        raw = read_ind(raw_path)
        filtered = read_ind(filtered_path)
        removed = sorted(set(raw) - set(filtered))
        added = sorted(set(filtered) - set(raw))
        evidence["panels"][panel_key] = {
            "raw_ind": str(raw_path),
            "filtered_ind": str(filtered_path),
            "raw_n": len(raw),
            "filtered_n": len(filtered),
            "removed_ids": removed,
            "added_ids": added,
            "raw_labels": sorted({row["label"] for row in raw.values()}),
            "filtered_labels": sorted({row["label"] for row in filtered.values()}),
        }
        for sample_id in removed:
            diff_rows.append({
                "panel": panel_key, "change": "removed_by_filtered",
                "sample_id": sample_id, "raw_label": raw[sample_id]["label"],
            })
        for sample_id in added:
            diff_rows.append({
                "panel": panel_key, "change": "added_by_filtered",
                "sample_id": sample_id, "raw_label": "NA",
            })

    bam_dir = Path(cfg["inputs"]["ancient_bam_dir"])
    for bam in sorted(bam_dir.glob("*.bam")):
        counts = samtools_flagstat(bam)
        evidence["ancient_bams"].append({
            "bam": str(bam), **counts,
        })

    with (out_dir / "server_evidence.json").open("w") as handle:
        json.dump(evidence, handle, indent=2, sort_keys=True)
        handle.write("\n")
    with (out_dir / "panel_filtered_sample_diff.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["panel", "change", "sample_id", "raw_label"], delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(diff_rows)
    with (out_dir / "ancient_bam_pair_flags.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["bam", "records", "paired_flag", "proper_pair_flag"], delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(evidence["ancient_bams"])

    print(f"wrote {out_dir / 'server_evidence.json'}")
    print(f"panel sample-diff rows: {len(diff_rows)}")
    print(f"ancient BAMs counted: {len(evidence['ancient_bams'])}")


if __name__ == "__main__":
    main()
