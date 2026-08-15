#!/usr/bin/env python3
"""
Aggregate run_masked_loo_validation.sh's per-individual
summarize_projection_distances.py outputs into a confusion matrix.

This does NOT implement any new classification rule -- "predicted
label" is simply whichever population summarize_projection_distances.py
already ranked #1 (nearest centroid) for that masked individual, exactly
the same ranking already used to describe real ancient samples. This
script only counts and cross-tabulates results that already exist on
disk.

Reads OUT_DIR/manifest.tsv (true_label, held_out_id, prefix) and each
row's "<prefix>.nearest.tsv" (rank-ordered population distances from
summarize_projection_distances.py), and reports:
  - a confusion matrix (true label vs. #1-ranked predicted label)
  - per-true-label accuracy (predicted #1 == true label)
  - overall accuracy
  - the aromatic<->japonica specific cross-confusion rate, if both
    labels are present in the manifest
  - for each individual, the rank at which its own true label appears
    in its ranking (not just whether #1 was correct) -- a near miss
    (true label ranked #2 close behind #1) is a different result than
    the true label being nowhere near the top
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--manifest", required=True, help="run_masked_loo_validation.sh's OUT_DIR/manifest.tsv")
    parser.add_argument("--out-matrix", default=None, help="optional TSV: confusion matrix")
    parser.add_argument("--out-detail", default=None, help="optional TSV: one row per individual with predicted label + rank of true label")
    return parser.parse_args(argv)


def load_nearest(path: Path) -> list[tuple[int, str, float]]:
    """Returns [(rank, population, distance), ...] sorted by rank."""
    rows = []
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows.append((int(row["rank"]), row["population"], float(row["distance"])))
    rows.sort(key=lambda r: r[0])
    return rows


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    manifest_path = Path(args.manifest)

    records = []
    with manifest_path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            true_label = row["true_label"]
            held_out_id = row["held_out_id"]
            prefix = Path(row["prefix"])
            nearest_path = prefix.parent / f"{prefix.name}.nearest.tsv"
            if not nearest_path.is_file():
                print(f"WARNING: {nearest_path} missing (from manifest row {true_label}/{held_out_id}) -- skipping, run may still be in progress", file=sys.stderr)
                continue
            ranking = load_nearest(nearest_path)
            if not ranking:
                print(f"WARNING: {nearest_path} has no ranked populations -- skipping", file=sys.stderr)
                continue
            predicted_label = ranking[0][1]
            predicted_dist = ranking[0][2]
            true_rank = next((r for r, lab, _ in ranking if lab == true_label), None)
            true_dist = next((d for _, lab, d in ranking if lab == true_label), None)
            records.append({
                "true_label": true_label,
                "held_out_id": held_out_id,
                "predicted_label": predicted_label,
                "predicted_distance": predicted_dist,
                "correct": predicted_label == true_label,
                "true_label_rank": true_rank,
                "true_label_distance": true_dist,
            })

    if not records:
        print("ERROR: no usable records found -- check manifest and .nearest.tsv paths", file=sys.stderr)
        return 1

    true_labels = sorted({r["true_label"] for r in records})
    predicted_labels = sorted({r["predicted_label"] for r in records})
    all_labels = sorted(set(true_labels) | set(predicted_labels))

    matrix: dict[str, dict[str, int]] = {t: defaultdict(int) for t in true_labels}
    for r in records:
        matrix[r["true_label"]][r["predicted_label"]] += 1

    print(f"\n[confusion matrix] rows = true label, columns = predicted (#1-ranked) label, n = {len(records)} individuals\n")
    header = "true\\pred".ljust(22) + "".join(lab[:12].rjust(14) for lab in all_labels) + "  total  accuracy"
    print(header)
    for t in true_labels:
        total = sum(matrix[t].values())
        correct = matrix[t].get(t, 0)
        acc = correct / total if total else float("nan")
        row_str = t.ljust(22) + "".join(str(matrix[t].get(lab, 0)).rjust(14) for lab in all_labels)
        row_str += f"  {total:5d}  {acc:6.1%}"
        print(row_str)

    total_correct = sum(1 for r in records if r["correct"])
    print(f"\n[overall accuracy] {total_correct}/{len(records)} = {total_correct / len(records):.1%}")

    if "aromatic" in true_labels and "japonica" in true_labels:
        aro_records = [r for r in records if r["true_label"] == "aromatic"]
        jap_records = [r for r in records if r["true_label"] == "japonica"]
        aro_as_jap = sum(1 for r in aro_records if r["predicted_label"] == "japonica")
        jap_as_aro = sum(1 for r in jap_records if r["predicted_label"] == "aromatic")
        print("\n[aromatic <-> japonica cross-confusion]")
        if aro_records:
            print(f"  true aromatic predicted japonica: {aro_as_jap}/{len(aro_records)} = {aro_as_jap / len(aro_records):.1%}")
        if jap_records:
            print(f"  true japonica predicted aromatic: {jap_as_aro}/{len(jap_records)} = {jap_as_aro / len(jap_records):.1%}")
        combined = aro_records + jap_records
        combined_confused = aro_as_jap + jap_as_aro
        if combined:
            print(f"  combined aromatic/japonica mutual confusion rate: {combined_confused}/{len(combined)} = {combined_confused / len(combined):.1%}")

    if args.out_matrix:
        out_path = Path(args.out_matrix)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with out_path.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["true_label"] + all_labels + ["total", "accuracy"])
            for t in true_labels:
                total = sum(matrix[t].values())
                correct = matrix[t].get(t, 0)
                acc = correct / total if total else float("nan")
                writer.writerow([t] + [matrix[t].get(lab, 0) for lab in all_labels] + [total, f"{acc:.4f}"])
        print(f"\n[done] wrote {out_path}", file=sys.stderr)

    if args.out_detail:
        out_path = Path(args.out_detail)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with out_path.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["true_label", "held_out_id", "predicted_label", "predicted_distance", "correct", "true_label_rank", "true_label_distance"])
            for r in records:
                writer.writerow([
                    r["true_label"], r["held_out_id"], r["predicted_label"],
                    f"{r['predicted_distance']:.6f}", r["correct"],
                    r["true_label_rank"] if r["true_label_rank"] is not None else "NA",
                    f"{r['true_label_distance']:.6f}" if r["true_label_distance"] is not None else "NA",
                ])
        print(f"[done] wrote {out_path}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
