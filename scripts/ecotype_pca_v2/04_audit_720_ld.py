#!/usr/bin/env python3
"""Phase 1 / Batch 1 (corrected), script 04.

LD decay by physical-distance bin, on a handful of chromosomes (default
1,3,6,12 per docs/ECOTYPE_PCA_PANEL_QC_DESIGN.md section 3 step 5). Written
for panel B (6.7M_720) since that panel's LD structure is the least
understood of the three, but usable generically via --bfile/--keep/--label.

Audit only -- this determines what LD-pruning window is appropriate, it
does not itself prune anything, and it does not change what is measured
(genome-wide pairwise r^2 within --max-window-kb, binned by physical
distance -- exactly what the original naive version computed).

CORRECTED (2026-08-15, GPT review of commit 10878d7): the original version
called `plink2 --r2-unphased --ld-window-kb 500 --ld-window-r2 0
--ld-window 999999` once per whole chromosome. On a dense panel like
6.7M_720 (~57 SNP/kb per docs/ECOTYPE_PCA_PANEL_QC_DESIGN.md section 3),
that is an uncontrolled, potentially enormous pairwise output -- both a
disk-space and a Python-memory risk (the old code appended every single
pair's r^2 into an in-memory list per distance bin). Neither problem
requires changing what is measured, so both are fixed here as pure
engineering/resource-safety issues (explicitly in scope per the batch
correction rules), not as new statistical parameters:

1. Genome-wide-in-one-shot -> processed in fixed-size physical CHUNKS
   (--block-mb, default 20Mb) with a halo of --max-window-kb behind each
   chunk's start, so every true pair within the window is still found
   (nothing is dropped, nothing is subsampled) but no single plink2
   invocation or intermediate .vcor file spans an entire chromosome. Each
   pair is attributed to exactly one chunk (the chunk owning its
   lower-position SNP) so nothing is double-counted. Chunk .vcor files are
   deleted after being parsed, bounding disk usage to ~1 chunk at a time.
2. Per-bin raw r^2 lists -> streaming (count, sum) accumulation. Reports
   the exact same mean r^2 per bin as before; just O(1) memory per bin
   instead of O(n_pairs).

--block-mb is an execution/chunking parameter only -- it changes how the
computation is scheduled, not what is computed (same window, same r2
threshold of 0 meaning "keep everything", same bins). It is NOT a
statistical parameter and was not invented to replace one; if chunking
alone still turns out to be operationally infeasible for the full 6.7M
panel, that is the point to stop and propose actual subsampling for
approval -- not something decided here.

Also corrected: any required chromosome that fails (plink2 error, missing
.vcor, incompatible header, or zero SNPs found in --bfile at all) now makes
the whole script exit non-zero. The previous version logged an error and
`continue`d, silently returning exit 0 with a partial/empty summary --
exactly the failure mode this correction batch was told to fix.
"""
import shutil
import subprocess
import sys
from bisect import bisect_left, bisect_right
from pathlib import Path
from collections import defaultdict

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_ecotype_v2 import base_argparser, setup_logger, load_config, check_output_not_present

DIST_BINS_BP = [
    (0, 1_000, "0-1kb"), (1_000, 2_000, "1-2kb"), (2_000, 5_000, "2-5kb"),
    (5_000, 10_000, "5-10kb"), (10_000, 20_000, "10-20kb"), (20_000, 50_000, "20-50kb"),
    (50_000, 100_000, "50-100kb"), (100_000, 200_000, "100-200kb"), (200_000, 500_000, "200-500kb"),
]


def bin_label(dist):
    for lo, hi, lab in DIST_BINS_BP:
        if lo <= dist < hi:
            return lab
    return None


def run_plink_ld(bfile, keep, chrom, from_bp, to_bp, max_window_kb, out_prefix, logger):
    cmd = ["plink2", "--bfile", bfile, "--chr", str(chrom),
           "--from-bp", str(from_bp), "--to-bp", str(to_bp),
           "--r2-unphased", "--ld-window-kb", str(max_window_kb),
           "--ld-window-r2", "0", "--ld-window", "999999",
           "--out", str(out_prefix)]
    if keep:
        cmd[3:3] = ["--keep", keep]
    logger.info("RUN: " + " ".join(cmd))
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.stderr.strip():
        logger.info("STDERR:\n" + proc.stderr.strip())
    return proc


def main():
    ap = base_argparser(__doc__)
    ap.add_argument("--label", required=True)
    ap.add_argument("--bfile", required=True)
    ap.add_argument("--keep", default=None)
    ap.add_argument("--chroms", default="1,3,6,12")
    ap.add_argument("--max-window-kb", type=int, default=500)
    ap.add_argument("--block-mb", type=float, default=20,
                     help="execution chunk size in Mb -- resource-safety knob only, "
                          "does not change what is measured (see script docstring)")
    ap.add_argument("--keep-chunk-files", action="store_true",
                     help="do not delete per-chunk plink2 output (debugging only)")
    args = ap.parse_args()
    logger, _ = setup_logger(f"04_audit_720_ld.{args.label}", args.out_dir)
    load_config(args.config)

    out_dir = Path(args.out_dir)
    summary_path = out_dir / f"{args.label}.ld_decay.summary.tsv"
    work_dir = out_dir / f"{args.label}.ld_decay_work"
    check_output_not_present([summary_path], args.overwrite, logger)
    if work_dir.exists() and any(work_dir.iterdir()) and not args.overwrite:
        logger.error(f"refusing to overwrite existing non-empty work dir (pass --overwrite): {work_dir}")
        sys.exit(3)
    if work_dir.exists() and args.overwrite:
        shutil.rmtree(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    chroms = [c.strip() for c in args.chroms.split(",") if c.strip()]
    block_bp = int(args.block_mb * 1_000_000)
    halo_bp = args.max_window_kb * 1000

    logger.info(f"chroms={chroms} block_bp={block_bp} halo_bp={halo_bp} "
                f"max_window_kb={args.max_window_kb}")

    positions_by_chrom = defaultdict(list)
    with open(args.bfile + ".bim") as fh:
        for line in fh:
            parts = line.split()
            chrom, pos = parts[0], int(parts[3])
            if chrom in chroms:
                positions_by_chrom[chrom].append(pos)
    for c in positions_by_chrom:
        positions_by_chrom[c].sort()

    binned = defaultdict(lambda: defaultdict(lambda: [0, 0.0]))  # chrom -> bin -> [n, sum]
    failed_chroms = []

    for chrom in chroms:
        positions = positions_by_chrom.get(chrom, [])
        if not positions:
            logger.error(f"chr{chrom}: FAIL -- zero SNPs found in {args.bfile}.bim for this "
                         f"chromosome. Required chromosome with no data is a hard failure, "
                         f"not silently skipped -- check --chroms / --bfile.")
            failed_chroms.append(chrom)
            continue

        chrom_max = positions[-1]
        n_blocks = chrom_max // block_bp + 1
        logger.info(f"chr{chrom}: {len(positions)} SNPs, max_pos={chrom_max}, "
                    f"{n_blocks} chunk(s) of {block_bp}bp")

        chrom_failed = False
        for block_i in range(n_blocks):
            block_start = block_i * block_bp + 1
            block_end = min(block_start + block_bp - 1, chrom_max)
            if bisect_right(positions, block_end) <= bisect_left(positions, block_start):
                continue  # no SNPs owned by this block, skip invoking plink2 entirely

            plink_from = max(1, block_start - halo_bp)
            chunk_prefix = work_dir / f"chr{chrom}.block{block_i}"
            proc = run_plink_ld(args.bfile, args.keep, chrom, plink_from, block_end,
                                 args.max_window_kb, chunk_prefix, logger)
            if proc.returncode != 0:
                logger.error(f"chr{chrom} block {block_i} [{block_start}-{block_end}]: FAIL -- "
                             f"plink2 exited {proc.returncode}. Likely a flag-syntax mismatch "
                             f"for this plink2 build (server has v2.0.0-a.6.9LM per prior "
                             f"sessions) -- do not substitute a different LD tool without "
                             f"reporting this first.")
                chrom_failed = True
                break

            vcor_path = None
            for cand in (chunk_prefix.with_suffix(".vcor"), Path(str(chunk_prefix) + ".vcor")):
                if cand.is_file():
                    vcor_path = cand
                    break
            if vcor_path is None:
                logger.error(f"chr{chrom} block {block_i}: FAIL -- expected .vcor output "
                             f"missing next to {chunk_prefix} (plink2 exited 0 but produced "
                             f"no output -- check {work_dir} for what it actually wrote)")
                chrom_failed = True
                break

            with open(vcor_path) as fh:
                header = fh.readline().lstrip("#").split()
                try:
                    pos_a_idx = header.index("POS_A")
                    pos_b_idx = header.index("POS_B")
                    r2_idx = header.index("UNPHASED_R2")
                except ValueError:
                    logger.error(f"chr{chrom} block {block_i}: FAIL -- incompatible .vcor "
                                 f"header {header} (expected POS_A/POS_B/UNPHASED_R2) -- "
                                 f"this plink2 build's output format differs from what this "
                                 f"script assumes, needs an implementation fix, not a re-run")
                    chrom_failed = True
                    break
                n_pairs_chunk = 0
                for line in fh:
                    parts = line.split()
                    if len(parts) <= max(pos_a_idx, pos_b_idx, r2_idx):
                        continue
                    pos_a, pos_b = int(parts[pos_a_idx]), int(parts[pos_b_idx])
                    owner_pos = min(pos_a, pos_b)
                    if not (block_start <= owner_pos <= block_end):
                        continue  # owned by a different (already- or not-yet-processed) block
                    dist = abs(pos_b - pos_a)
                    lab = bin_label(dist)
                    if lab is None:
                        continue
                    entry = binned[chrom][lab]
                    entry[0] += 1
                    entry[1] += float(parts[r2_idx])
                    n_pairs_chunk += 1
            logger.info(f"chr{chrom} block {block_i} [{block_start}-{block_end}]: "
                        f"{n_pairs_chunk} owned pairs")

            if not args.keep_chunk_files:
                for p in work_dir.glob(f"chr{chrom}.block{block_i}.*"):
                    p.unlink()

        if chrom_failed:
            failed_chroms.append(chrom)

    rows = []
    for chrom in chroms:
        for lo, hi, lab in DIST_BINS_BP:
            n, s = binned[chrom].get(lab, [0, 0.0])
            mean_r2 = (s / n) if n else None
            rows.append((chrom, lab, n, mean_r2))

    with open(summary_path, "w") as fh:
        fh.write("chrom\tdistance_bin\tn_pairs\tmean_r2\n")
        for chrom, lab, n, mean_r2 in rows:
            fh.write(f"{chrom}\t{lab}\t{n}\t{'' if mean_r2 is None else f'{mean_r2:.5f}'}\n")
    logger.info(f"wrote {summary_path} (diagnostic output written regardless of failure below)")
    for chrom, lab, n, mean_r2 in rows:
        logger.info(f"chr{chrom} {lab}: n={n} mean_r2={mean_r2}")

    if not args.keep_chunk_files:
        try:
            work_dir.rmdir()
        except OSError:
            pass  # non-empty (some debug leftovers) -- harmless, not worth failing over

    if failed_chroms:
        logger.error(f"FAIL: required chromosome(s) failed: {failed_chroms}. "
                     f"summary.tsv above is partial/diagnostic only, NOT a passing result.")
        sys.exit(2)
    sys.exit(0)


if __name__ == "__main__":
    main()
