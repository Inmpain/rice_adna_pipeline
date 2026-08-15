#!/usr/bin/env python3
"""
Plot PC1 vs PC2 from one or more smartpca .evec files, colored by
population label, with named samples (ancient/simulated) drawn as
highlighted stars.

WARNING -- axes are NOT comparable across different .evec files unless
they were built with the exact same poplistname (same reference
individuals defining the axes). Two runs whose axis-building reference
set differs (e.g. "wild rice included" vs "wild rice excluded") produce
independently rotated/scaled coordinate systems: PC1/PC2 values and
inter-point distances from one run cannot be compared numerically
against another run, only qualitative structure can (which population a
sample lands nearest to, relative cluster separation). This is exactly
why each --evec gets its own subplot with its own independent axis
scale, never a shared/overlaid plot. This applies doubly to
run_sample_panel_pca.sh batch output (ECOTYPE_PCA_PHASE1_COMMANDS.md
section 7): each ancient sample there also gets its own private marker
SUBSET (shrunk to that sample's own covered SNPs), on top of whatever
reference set built the axes -- so even two civan-panel subplots in the
same grid are two independently-built PCAs, not just two views of one
shared space. Never conclude "sample A's PC1 is bigger than sample B's
PC1" from this plot; only conclude "sample A projects closer to
population X than to population Y" within its own subplot.

Usage (single panel):
  python3 plot_pca_projection.py \\
    --evec civan_refonly=LV7008416379.civan.TV.evec \\
    --highlight LV7008416379 \\
    --out civan_refonly.png

Usage (side-by-side before/after comparison):
  python3 plot_pca_projection.py \\
    --evec wild_in_axis=loo_smoke/civan_loo_test.evec \\
    --evec domesticated_only=civan_refonly_check/LV7008416379.civan.TV.evec \\
    --highlight LV7008416379 \\
    --title "Civan panel: wild-rice-in-axis fix, before vs after" \\
    --out civan_before_after.png

Usage (batch grid -- e.g. all 16 samples for one panel from
run_sample_panel_pca.sh's scale-out, ECOTYPE_PCA_PHASE1_COMMANDS.md
section 7): each matched file becomes its own subplot, labeled by the
part of its filename before the first '.' (i.e. the sample ID for
run_sample_panel_pca.sh's SAMPLE.PANEL.TRACK.evec naming). No
--highlight needed for these -- the ancient sample in each file already
carries the "Ancient" population label and is auto-highlighted.
  python3 plot_pca_projection.py \\
    --evec-glob "/path/to/pca_runs/*.civan.TV.evec" \\
    --ncols 4 \\
    --title "Civan panel, all 16 ancient samples (each subplot its own axis scale)" \\
    --out civan_all16_grid.png
"""

from __future__ import annotations

import argparse
import glob
import sys
from collections import defaultdict
from pathlib import Path


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--evec",
        action="append",
        default=[],
        help="PANEL_LABEL=path/to/file.evec (repeatable) -- one subplot per entry, in the order given",
    )
    parser.add_argument(
        "--evec-glob",
        action="append",
        default=[],
        help="glob pattern (repeatable), e.g. 'pca_runs/*.civan.TV.evec' -- each match becomes its own "
        "subplot, in sorted filename order, labeled by the filename's first '.'-delimited field "
        "(the sample ID for run_sample_panel_pca.sh output); added after any --evec entries",
    )
    parser.add_argument(
        "--highlight",
        action="append",
        default=[],
        help="individual ID to draw as a highlighted star and annotate by name (repeatable) -- "
        "usually unnecessary since any row labeled 'Ancient' is already auto-highlighted",
    )
    parser.add_argument("--pc-x", type=int, default=1, help="which PC (1-indexed) to plot on the x axis (default 1)")
    parser.add_argument("--pc-y", type=int, default=2, help="which PC (1-indexed) to plot on the y axis (default 2)")
    parser.add_argument("--ncols", type=int, default=4, help="subplot grid columns (default 4); rows added as needed")
    parser.add_argument("--min-pop-size", type=int, default=5, help="labels with fewer individuals (summed across all given evec files) are pooled into a single gray 'other (n<N)' bucket")
    parser.add_argument("--title", default=None, help="overall figure title")
    parser.add_argument("--out", required=True, help="output image path (.png)")
    args = parser.parse_args(argv)
    if not args.evec and not args.evec_glob:
        parser.error("at least one --evec or --evec-glob is required")
    return args


def load_evec(path: Path, pc_x: int, pc_y: int) -> list[tuple[str, float, float, str]]:
    """pc_x/pc_y are 1-indexed PC numbers (PC1, PC2, ...), matching .evec's own column numbering."""
    rows = []
    with path.open() as handle:
        first = handle.readline()
        if not first.strip().startswith("#"):
            raise ValueError(f"{path}: expected first line to start with '#' (smartpca's eigvals header, possibly indented), got {first[:40]!r}")
        for line_no, line in enumerate(handle, start=2):
            fields = line.split()
            need = 1 + max(pc_x, pc_y) + 1  # id + PCs up to the higher requested one + label
            if len(fields) < need:
                raise ValueError(f"{path}:{line_no}: expected at least {need} fields (id, PC1..PC{max(pc_x, pc_y)}, ..., label) to plot PC{pc_x}/PC{pc_y}, got {len(fields)}: {line!r}")
            sample_id = fields[0]
            x = float(fields[pc_x])
            y = float(fields[pc_y])
            label = fields[-1]
            rows.append((sample_id, x, y, label))
    return rows


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    panels: list[tuple[str, list[tuple[str, float, float, str]]]] = []
    for spec in args.evec:
        if "=" in spec:
            panel_label, path_str = spec.split("=", 1)
        else:
            panel_label, path_str = Path(spec).stem, spec
        panels.append((panel_label, load_evec(Path(path_str), args.pc_x, args.pc_y)))

    for pattern in args.evec_glob:
        matches = sorted(glob.glob(pattern))
        if not matches:
            print(f"ERROR: --evec-glob pattern matched no files: {pattern!r}", file=sys.stderr)
            return 1
        for path_str in matches:
            panel_label = Path(path_str).name.split(".")[0]
            panels.append((panel_label, load_evec(Path(path_str), args.pc_x, args.pc_y)))

    counts: dict[str, int] = defaultdict(int)
    for _, rows in panels:
        for _, _, _, label in rows:
            counts[label] += 1

    kept_labels = sorted(
        label for label, n in counts.items()
        if n >= args.min_pop_size and label not in ("Ancient", "LOO_HELDOUT_EXCLUDED")
    )
    cmap = plt.get_cmap("tab20")
    color_of = {label: cmap(i % 20) for i, label in enumerate(kept_labels)}
    other_color = (0.55, 0.55, 0.55, 0.45)

    ncols = max(1, min(args.ncols, len(panels)))
    nrows = -(-len(panels) // ncols)  # ceil division
    fig, axes = plt.subplots(nrows, ncols, figsize=(6.2 * ncols, 5.6 * nrows), squeeze=False)
    axes_flat = axes.flatten()

    highlight_set = set(args.highlight)

    for ax, (panel_label, rows) in zip(axes_flat, panels):
        legend_seen: set[str] = set()
        for sample_id, pc1, pc2, label in rows:
            if sample_id in highlight_set or label in ("Ancient", "LOO_HELDOUT_EXCLUDED"):
                continue
            color = color_of.get(label, other_color)
            plot_label = None
            if label in color_of and label not in legend_seen:
                plot_label = label
                legend_seen.add(label)
            elif label not in color_of and "other (n<{})".format(args.min_pop_size) not in legend_seen:
                plot_label = f"other (n<{args.min_pop_size})"
                legend_seen.add(plot_label)
            ax.scatter(pc1, pc2, s=14, color=color, alpha=0.65, linewidths=0, label=plot_label)

        for sample_id, pc1, pc2, label in rows:
            if sample_id not in highlight_set and label not in ("Ancient", "LOO_HELDOUT_EXCLUDED"):
                continue
            ax.scatter(pc1, pc2, s=170, marker="*", color="red", edgecolor="black", linewidths=0.9, zorder=5)
            ax.annotate(
                f"{sample_id} ({label})" if label not in ("Ancient",) else sample_id,
                (pc1, pc2), textcoords="offset points", xytext=(7, 7),
                fontsize=8, fontweight="bold",
            )

        ax.set_xlabel(f"PC{args.pc_x} (this panel's own axis scale)")
        ax.set_ylabel(f"PC{args.pc_y} (this panel's own axis scale)")
        ax.set_title(panel_label)
        ax.legend(fontsize=6, markerscale=1.4, loc="best", ncol=1, framealpha=0.85)

    for ax in axes_flat[len(panels):]:
        ax.axis("off")

    if args.title:
        fig.suptitle(args.title)
    fig.tight_layout()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    print(f"[done] wrote {out_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
