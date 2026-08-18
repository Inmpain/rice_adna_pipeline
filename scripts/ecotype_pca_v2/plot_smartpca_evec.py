#!/usr/bin/env python3
"""Plot a smartpca .evec/.eval pair: modern samples colored by population label,
ancient samples (label 'Ancient') as normal black points with name labels (no star).
Axis labels show per-PC explained variance; title shows marker count."""
import argparse
from collections import defaultdict


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--evec", required=True)
    ap.add_argument("--eval", required=True)
    ap.add_argument("--ind", required=True)
    ap.add_argument("--nmarkers", type=int, required=True)
    ap.add_argument("--title", default="")
    ap.add_argument("--out-prefix", required=True)
    return ap.parse_args()


def load_eval(path):
    vals = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line:
                vals.append(float(line))
    return vals


def load_evec(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.lstrip().startswith("#"):
                continue
            f = line.split()
            if len(f) < 3:
                continue
            iid = f[0]
            label = f[-1]
            vals = list(map(float, f[1:-1]))
            rows.append((iid, vals, label))
    return rows


def main():
    args = parse_args()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    evals = load_eval(args.eval)
    total = sum(evals) or 1.0
    rows = load_evec(args.evec)
    npc = len(rows[0][1])
    by_lab = defaultdict(list)
    ancient = []
    for iid, vals, label in rows:
        if label == "Ancient":
            ancient.append((iid, vals))
        else:
            by_lab[label].append((iid, vals))

    npairs = min(5, npc // 2)
    fig, axes = plt.subplots(1, npairs, figsize=(6.8 * npairs, 5.8), squeeze=False)
    cmap = plt.get_cmap("tab20")
    labs = sorted(by_lab)
    for pi in range(npairs):
        ax = axes[0][pi]
        xi, yi = 2 * pi, 2 * pi + 1
        for li, lab in enumerate(labs):
            xs = [v[xi] for _, v in by_lab[lab]]
            ys = [v[yi] for _, v in by_lab[lab]]
            ax.scatter(xs, ys, s=10, color=cmap(li % 20), label=lab, alpha=0.6, linewidths=0)
        for iid, v in ancient:
            ax.scatter(v[xi], v[yi], s=70, color="black", marker="o",
                       edgecolor="white", linewidths=0.5, zorder=5)
            ax.annotate(iid, (v[xi], v[yi]), textcoords="offset points",
                        xytext=(5, 5), fontsize=7, fontweight="bold")
        ax.set_xlabel(f"PC{xi + 1} ({evals[xi] / total * 100:.2f}%)")
        ax.set_ylabel(f"PC{yi + 1} ({evals[yi] / total * 100:.2f}%)")
        ax.set_title(f"{args.title} | markers={args.nmarkers}")
        ax.legend(fontsize=5, markerscale=1.5, loc="best")
    fig.tight_layout()
    out = f"{args.out_prefix}.png"
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
