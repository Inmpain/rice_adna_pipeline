#!/usr/bin/env python3
"""
Per-read Oryza vs non-Oryza competitive best-hit filter, with ancient-DNA
terminal-damage-aware NM correction.

v2 (2026-08-08) -- see docs/ORYZA_BESTHIT_HANDOFF.md section 5.1b for the full
design writeup. Two changes vs v1 (archived as
oryza_besthit_damage_filter_v1.py, kept for reference/rollback):

1. Oryza scope is now the WHOLE genus by default: every species-rank
   descendant of taxid 4527 (Oryza) in nodes.dmp, resolved dynamically at
   startup (see Taxonomy.genus_species_taxids below) -- not a hardcoded
   3-species list (rufipogon/sativa/nivara). This means a read that is a
   genuine match to, say, O. longistaminata (well represented in our WGS
   shard, ~905 contigs) is no longer forced to compete against sativa/
   rufipogon/nivara as if it were foreign contamination; it's Oryza either
   way. --oryza-taxids still exists as an explicit manual override (pass it
   to reproduce v1's narrower 3-species behavior exactly); when omitted,
   --oryza-genus-taxid (default 4527) drives the auto-resolution.

2. Optional (OFF BY DEFAULT) quality pre-gate, --min-best-similarity /
   --max-best-raw-nm: reject a read outright, before any per-species
   classification, if its single best raw-NM hit across ALL species fails
   both thresholds (OR-gated if both given, matching the source idea).
   Borrowed from a labmate's classifier
   (besthit_competitive_top10_showOryza_optimized.py's --target-min-sim /
   --target-max-nm). Off by default because it changes which reads even
   reach the Oryza-vs-non-Oryza competition and has not been validated
   against our data.

One thing considered and DELIBERATELY NOT adopted from that labmate script:
its faster MD-tag-string damage parsing also skips damage evaluation for any
alignment whose RAW NM isn't already that species' minimum -- which assumes
raw NM ranks the same as damage-adjusted NM. That is not guaranteed: an
alignment with a *higher* raw NM can still end up with a *lower* adjusted NM
if more of its mismatches happen to land inside the terminal-damage window.
That risk is bigger for us than for the labmate script, because our default
--damage-window is 5bp per end (up to 10 possible credit) vs their 1bp (up to
2). alignment_metrics() below is therefore left exactly as validated in v1
(full pysam.get_aligned_pairs(with_seq=True) walk over EVERY alignment, no
raw-NM pre-pruning) rather than risk a silent regression from an unvalidated
rewrite. Revisit if per-sample runtime actually becomes a bottleneck -- with
a proper side-by-side check against v1's output on real data first.

======================================================================
Per-alignment metrics (needs the BAM's NM tag, MD tag and query SEQ):
======================================================================
  NM                  : bowtie2's own NM tag (edit distance: subs + indels).
  substitution_count   : mismatches only, read off MD/CIGAR (via
                         pysam.get_aligned_pairs(with_seq=True); a lowercase
                         ref base marks a mismatch). Indels excluded.
  terminal_damage_count: substitutions within --damage-window bases of a read
                         end that match the aDNA deamination signature, IN THE
                         READ'S OWN 5'->3' ORIENTATION:
                             5' end: ref=C, read=T
                             3' end: ref=G, read=A
                         For a reverse-strand alignment the BAM's SEQ is
                         already the reverse complement of the sequenced read
                         (SAM convention), so recovering "the read's own
                         orientation" means reversing the aligned (ref,read)
                         base pairs AND complementing both bases -- a plain
                         reverse (no complement) would leave the reference-
                         strand G->A signature in place instead of turning it
                         back into the read-strand C->T signature the rule
                         above expects.
  adjusted_NM           : NM - terminal_damage_count (indels are never
                         subtracted; only genuine terminal deamination subs).

If an alignment has no usable SEQ/MD (query_sequence is None, or MD is
absent), substitution_count/terminal_damage_count cannot be computed for that
alignment; it falls back to adjusted_NM = NM (no damage credit) and is
counted in the [warn] reads_missing_md_or_seq stat printed at the end -- this
should be ~0 for a normal bowtie2 BAM; investigate if it isn't.

======================================================================
Per-read decision:
======================================================================
Alignments are grouped by SPECIES taxid (BAM reference -> acc2taxid -> walk
up nodes.dmp to the nearest rank=="species" ancestor). Within a species, only
the single best alignment survives, ranked by
    (adjusted_NM asc, NM asc, substitution_count asc, AS desc, reference_name)
Non-Oryza species are ranked the same way and the best 10 are kept
(--top-n); the Oryza species (genus-wide by default, see above) are ALWAYS
kept in addition, never competing for one of the 10 non-Oryza slots.

  best_nonoryza = best-ranked non-Oryza species for this read (or none)
  best_oryza    = best-ranked Oryza species for this read (or none)

  no Oryza hit at all                                -> REJECT (no_oryza_hit)
  Oryza hit, no non-Oryza hit at all                  -> KEEP  (oryza_only_no_competitor)
  best_oryza.adjusted_NM <= best_nonoryza.adjusted_NM -> KEEP  (oryza_at_least_as_good)
  best_oryza.adjusted_NM >  best_nonoryza.adjusted_NM -> REJECT (nonoryza_better)

A read with alignments in the BAM but where EVERY alignment's reference
fails to resolve to a species taxid gets decision=UNCLASSIFIED (still gets a
decisions.tsv row, with reason=no_resolvable_species, for debuggability) and
is NOT written to the top10 audit table (nothing to rank). A candidate read
that never appears in the BAM at all (bowtie2 --no-unal dropped it -- no hit
anywhere in the 131 databases) gets no decisions.tsv row but IS folded into
the unclassified_reads bucket of the sample summary, so that
    input_reads == kept_reads + rejected_nonoryza_better + rejected_no_oryza
                   + rejected_low_quality + unclassified_reads
holds exactly (checked at the end, hard failure unless --limit-reads is set).
rejected_low_quality is always 0 unless the optional pre-gate (see above) is
enabled.

MAPQ is never used -- under bowtie2 -k 100 competitive mapping MAPQ is not a
meaningful best-hit signal (see docs/ORYZA_BESTHIT_HANDOFF.md).

======================================================================
Outputs (all written to a .tmp path first, then os.rename'd -- atomic):
======================================================================
  <outdir>/<sample>.besthit.top10_species.tsv.gz
  <outdir>/<sample>.oryza_filter.decisions.tsv.gz
  <outdir>/<sample>.besthit_oryza.fastq.gz
  <outdir>/<sample>.summary.tsv          (one-sample row; concatenate across
                                           samples for besthit_summary.tsv --
                                           done by the submit script, not here,
                                           so parallel per-sample SLURM jobs
                                           never race on one shared file.
                                           NOTE: this v2 summary.tsv has one
                                           more column, rejected_low_quality,
                                           than v1's -- don't `merge` v1 and
                                           v2 sample summaries together
                                           without checking columns line up.)
  <outdir>/<sample>.finished             (touched only after every output
                                           above succeeded AND the consistency
                                           check passed)
"""

import argparse
import gzip
import os
import shutil
import subprocess
import sys

import pysam

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


# ===========================================================================
# Taxonomy: nodes.dmp (parent + rank) / names.dmp (scientific name)
# ===========================================================================
class Taxonomy:
    def __init__(self, nodes_file, names_file):
        self.parent = {}
        self.rank = {}
        self.name = {}
        self.children = {}       # NEW in v2: parent_taxid -> set(child_taxid)
        self._species_cache = {}
        with open(nodes_file) as fh:
            for line in fh:
                p = line.split("|")
                if len(p) < 3:
                    continue
                tid = p[0].strip()
                par = p[1].strip()
                self.parent[tid] = par
                self.rank[tid] = p[2].strip()
        with open(names_file) as fh:
            for line in fh:
                if "scientific name" not in line:
                    continue
                p = line.split("|")
                if len(p) < 2:
                    continue
                self.name[p[0].strip()] = p[1].strip()
        for tid, par in self.parent.items():
            if tid == par:      # root points to itself; not its own child
                continue
            self.children.setdefault(par, set()).add(tid)

    def species_of(self, taxid):
        """Walk up from taxid to the nearest rank=='species' ancestor
        (or taxid itself). None if no species-rank ancestor exists before
        the root (e.g. taxid is only resolved to genus or higher)."""
        cached = self._species_cache.get(taxid)
        if cached is not None:
            return cached if cached != "" else None
        cur = taxid
        seen = set()
        result = None
        while cur in self.parent and cur not in seen:
            seen.add(cur)
            if self.rank.get(cur) == "species":
                result = cur
                break
            par = self.parent[cur]
            if par == cur:
                break
            cur = par
        for t in seen:
            self._species_cache[t] = result if result is not None else ""
        return result

    def sci_name(self, taxid):
        return self.name.get(taxid, taxid)

    def genus_species_taxids(self, genus_taxid):
        """NEW in v2. All species-rank taxids in the subtree rooted at
        genus_taxid -- a BFS/DFS over the children map, so it only visits
        this genus's own descendants (tens of nodes for Oryza), not the
        whole taxonomy dump (millions of nodes)."""
        if genus_taxid not in self.parent:
            return set()
        seen = set()
        stack = [genus_taxid]
        species = set()
        while stack:
            cur = stack.pop()
            if cur in seen:
                continue
            seen.add(cur)
            if self.rank.get(cur) == "species":
                species.add(cur)
            stack.extend(self.children.get(cur, ()))
        return species


def load_acc2taxid(path):
    print(f"[load] acc2taxid: {path}", file=sys.stderr)
    acc_map = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("accession"):
                continue
            parts = line.strip().split()
            if len(parts) >= 3:
                tax_id = parts[2]
                acc_map[parts[0]] = tax_id
                if parts[1] != parts[0]:
                    acc_map[parts[1]] = tax_id
    return acc_map


def build_refid2species(bam, acc2taxid, tax):
    references = bam.references
    refid2species = []
    n_no_taxid = n_no_species = 0
    for name in references:
        tid = acc2taxid.get(name) or acc2taxid.get(name.split(".")[0])
        if tid is None:
            n_no_taxid += 1
            refid2species.append(None)
            continue
        sp = tax.species_of(tid)
        if sp is None:
            n_no_species += 1
        refid2species.append(sp)
    print(f"[load] {len(references)} reference contigs: "
          f"{n_no_taxid} without any acc2taxid hit, "
          f"{n_no_species} resolved to a taxid but no species-rank ancestor",
          file=sys.stderr)
    return refid2species


# ===========================================================================
# Per-alignment metrics (unchanged from v1 -- see module docstring for why)
# ===========================================================================
def alignment_metrics(aln, window):
    """-> (NM, substitution_count_or_None, terminal_damage_count, AS)"""
    try:
        nm = aln.get_tag("NM")
    except KeyError:
        nm = 0
    try:
        as_score = aln.get_tag("AS")
    except KeyError:
        as_score = -10**9

    if aln.query_sequence is None:
        return nm, None, 0, as_score

    try:
        pairs = aln.get_aligned_pairs(matches_only=False, with_seq=True)
    except ValueError:
        # no MD tag on this alignment
        return nm, None, 0, as_score

    seq = aln.query_sequence
    ordered = []          # (ref_base_upper, read_base_upper, is_mismatch), ref-forward order
    for qpos, rpos, ref_base in pairs:
        if qpos is None or rpos is None:
            continue       # pure insertion or deletion: no (ref,read) base pair
        ordered.append((ref_base.upper(), seq[qpos].upper(), ref_base.islower()))

    substitution_count = sum(1 for _, _, mm in ordered if mm)

    if aln.is_reverse:
        # BAM SEQ is revcomp(original read); reverse+complement to recover
        # the read's own 5'->3' orientation (see module docstring).
        read_oriented = [(b.translate(COMPLEMENT), r.translate(COMPLEMENT))
                         for b, r, _ in reversed(ordered)]
    else:
        read_oriented = [(b, r) for b, r, _ in ordered]

    n = len(read_oriented)
    damage = 0
    for i, (ref_b, read_b) in enumerate(read_oriented):
        if i < window and ref_b == "C" and read_b == "T":
            damage += 1
        if i >= n - window and ref_b == "G" and read_b == "A":
            damage += 1

    return nm, substitution_count, damage, as_score


def rank_key(rec):
    """Ascending = better. rec = (adjusted_NM, NM, substitution_count, AS, ref_name)."""
    adj_nm, nm, sub, as_score, ref_name = rec
    sub_key = sub if sub is not None else nm  # unknown subs: don't let it look artificially good
    return (adj_nm, nm, sub_key, -as_score, ref_name)


def raw_nm_tag(aln):
    """NEW in v2 -- cheap tag-only NM lookup for the optional quality pre-gate,
    deliberately not routed through alignment_metrics() (no MD/CIGAR parsing
    needed just to find the single best raw NM across all of a read's hits)."""
    try:
        return aln.get_tag("NM")
    except KeyError:
        return 0


# ===========================================================================
# FASTQ streaming (candidate reads in, KEEP reads out)
# ===========================================================================
def count_fastq_reads(path):
    n = 0
    opener = (lambda p: gzip.open(p, "rt")) if path.endswith(".gz") else open
    with opener(path) as fh:
        for _ in fh:
            n += 1
    return n // 4


def extract_fastq(src_fastq, keep_names, out_path, threads):
    tmp_path = out_path + ".tmp"
    gz = shutil.which("pigz") or shutil.which("gzip")
    out_args = [gz, "-c"] + (["-p", str(threads)] if "pigz" in gz and threads > 1 else [])
    fout = open(tmp_path, "wb")
    out_proc = subprocess.Popen(out_args, stdin=subprocess.PIPE, stdout=fout)

    in_proc = None
    if src_fastq.endswith(".gz"):
        in_proc = subprocess.Popen([gz, "-dc", src_fastq], stdout=subprocess.PIPE,
                                   bufsize=1 << 22)
        stream = in_proc.stdout
    else:
        stream = open(src_fastq, "rb")

    written = set()
    n_out = 0
    w = out_proc.stdin.write
    try:
        while True:
            h = stream.readline()
            if not h:
                break
            seq = stream.readline()
            plus = stream.readline()
            qual = stream.readline()
            if not qual:
                break
            rid = h[1:].split(None, 1)[0].decode("ascii", "replace")
            if rid in keep_names and rid not in written:
                w(h); w(seq); w(plus); w(qual)
                written.add(rid)
                n_out += 1
    finally:
        stream.close()
        if in_proc:
            in_proc.wait()
        out_proc.stdin.close()
        out_proc.wait()
        fout.close()

    os.rename(tmp_path, out_path)
    return n_out, written


# ===========================================================================
# Main
# ===========================================================================
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--bam", required=True, help="query-name sorted competitive BAM")
    ap.add_argument("--fastq", required=True, help="original candidate FASTQ(.gz)")
    ap.add_argument("--acc2taxid", required=True)
    ap.add_argument("--nodes", required=True)
    ap.add_argument("--names", required=True)
    ap.add_argument("--oryza-taxids", nargs="+", default=None,
                    help="Explicit whitelist override. If given, used verbatim "
                         "instead of genus-wide auto-resolution -- pass "
                         "'4529 4530 4536' to reproduce v1's rufipogon/sativa/"
                         "nivara-only behavior. Omit (default) for genus-wide.")
    ap.add_argument("--oryza-genus-taxid", default="4527",
                    help="Used only when --oryza-taxids is not given: every "
                         "species-rank descendant of this taxid in nodes.dmp "
                         "becomes the Oryza whitelist. Default 4527 = genus "
                         "Oryza (all ~18+ species, not just the 3 focal "
                         "ones -- see module docstring).")
    ap.add_argument("--damage-window", type=int, default=5)
    ap.add_argument("--top-n", type=int, default=10,
                    help="non-Oryza species kept in the audit table per read")
    ap.add_argument("--min-best-similarity", type=float, default=None,
                    help="Optional pre-gate, OFF BY DEFAULT: reject a read "
                         "outright (before any per-species classification) if "
                         "its single best raw-NM hit across ALL species fails "
                         "this AND (if --max-best-raw-nm is also given) that "
                         "threshold too -- i.e. OR-gated, pass if either "
                         "condition holds. Borrowed idea from a labmate's "
                         "classifier's --target-min-sim; unvalidated against "
                         "our data, see docs/ORYZA_BESTHIT_HANDOFF.md.")
    ap.add_argument("--max-best-raw-nm", type=int, default=None,
                    help="See --min-best-similarity.")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--threads", type=int, default=2, help="BAM/FASTQ IO threads")
    ap.add_argument("--limit-reads", type=int, default=None,
                    help="smoke test: stop after this many distinct read names "
                         "from the BAM; relaxes the final consistency check")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    smoke = args.limit_reads is not None
    if smoke:
        print(f"[smoke-test] limiting to first {args.limit_reads} reads from the BAM",
              file=sys.stderr)

    tax = Taxonomy(args.nodes, args.names)
    acc2taxid = load_acc2taxid(args.acc2taxid)

    if args.oryza_taxids:
        for t in args.oryza_taxids:
            if t not in tax.parent:
                sys.exit(f"[error] --oryza-taxids {t} not found in nodes.dmp")
        oryza_set = set(args.oryza_taxids)
        print(f"[config] Oryza scope: manual override, {len(oryza_set)} taxid(s): "
              + ", ".join(f"{t}({tax.sci_name(t)})" for t in sorted(oryza_set)),
              file=sys.stderr)
    else:
        if args.oryza_genus_taxid not in tax.parent:
            sys.exit(f"[error] --oryza-genus-taxid {args.oryza_genus_taxid} "
                     f"not found in nodes.dmp")
        genus_rank = tax.rank.get(args.oryza_genus_taxid)
        if genus_rank != "genus":
            print(f"[warn] taxid {args.oryza_genus_taxid} has rank "
                  f"'{genus_rank}', not 'genus' -- proceeding anyway, but "
                  f"double check this is the taxid you meant", file=sys.stderr)
        oryza_set = tax.genus_species_taxids(args.oryza_genus_taxid)
        if not oryza_set:
            sys.exit(f"[error] no species-rank descendants found under taxid "
                     f"{args.oryza_genus_taxid} in nodes.dmp -- check the "
                     f"taxid or pass --oryza-taxids explicitly")
        print(f"[config] Oryza scope: whole genus (taxid {args.oryza_genus_taxid}, "
              f"{tax.sci_name(args.oryza_genus_taxid)}), {len(oryza_set)} species "
              f"resolved: "
              + ", ".join(f"{t}({tax.sci_name(t)})" for t in sorted(oryza_set)),
              file=sys.stderr)

    quality_gate_enabled = (args.min_best_similarity is not None or
                            args.max_best_raw_nm is not None)
    if quality_gate_enabled:
        print(f"[config] quality pre-gate ON: min_best_similarity="
              f"{args.min_best_similarity} max_best_raw_nm={args.max_best_raw_nm} "
              f"(OR-gated if both given)", file=sys.stderr)

    bam = pysam.AlignmentFile(args.bam, "rb", threads=max(1, args.threads))
    hd = bam.header.to_dict().get("HD", {})
    if isinstance(hd, list):
        hd = hd[0]
    if hd.get("SO") != "queryname":
        sys.exit(f"[error] BAM sort order is '{hd.get('SO')}', need queryname "
                 f"(samtools sort -n)")
    refid2species = build_refid2species(bam, acc2taxid, tax)
    del acc2taxid

    top10_tmp = os.path.join(args.outdir, f"{args.sample}.besthit.top10_species.tsv.gz.tmp")
    dec_tmp = os.path.join(args.outdir, f"{args.sample}.oryza_filter.decisions.tsv.gz.tmp")
    top10_fh = gzip.open(top10_tmp, "wt")
    dec_fh = gzip.open(dec_tmp, "wt")
    top10_fh.write("read_name\tread_length\tspecies_hit_count\ttop10_rank\t"
                   "is_always_included\tspecies_taxid\tspecies_name\treference_name\t"
                   "NM\tsubstitution_count\tterminal_damage_count\tadjusted_NM\tAS\n")
    dec_fh.write("read_name\tbest_nonoryza_taxid\tbest_nonoryza_name\t"
                "best_nonoryza_NM\tbest_nonoryza_damage\tbest_nonoryza_adjusted_NM\t"
                "best_oryza_taxid\tbest_oryza_name\tbest_oryza_NM\tbest_oryza_damage\t"
                "best_oryza_adjusted_NM\tdecision\treason\n")

    keep_names = set()
    n_reads_with_alignment = 0
    n_unclassified_in_bam = 0
    n_kept = 0
    n_rejected_nonoryza_better = 0
    n_rejected_no_oryza = 0
    n_rejected_low_quality = 0
    n_missing_md_or_seq = 0

    def flush_read(qname, qlen, alns):
        """alns: list of pysam.AlignedSegment for this one read."""
        nonlocal n_reads_with_alignment, n_unclassified_in_bam, n_kept
        nonlocal n_rejected_nonoryza_better, n_rejected_no_oryza, n_missing_md_or_seq
        nonlocal n_rejected_low_quality

        n_reads_with_alignment += 1

        if quality_gate_enabled:
            best_raw_nm = min(raw_nm_tag(a) for a in alns)
            best_sim = (100.0 * (qlen - best_raw_nm) / qlen) if qlen > 0 else 0.0
            checks = []
            if args.min_best_similarity is not None:
                checks.append(best_sim >= args.min_best_similarity)
            if args.max_best_raw_nm is not None:
                checks.append(best_raw_nm <= args.max_best_raw_nm)
            if not any(checks):
                n_rejected_low_quality += 1
                dec_fh.write(f"{qname}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t"
                            f"REJECT\tlow_quality_pregate\n")
                return

        best_by_species = {}   # species_taxid -> (rank_key_tuple, ref_name, NM, sub, dmg, adj, AS)
        for aln in alns:
            sp = refid2species[aln.reference_id]
            if sp is None:
                continue
            nm, sub, dmg, as_score = alignment_metrics(aln, args.damage_window)
            if sub is None:
                n_missing_md_or_seq += 1
            adj = nm - dmg
            ref_name = bam.get_reference_name(aln.reference_id)
            rec = (adj, nm, sub, as_score, ref_name)
            prev = best_by_species.get(sp)
            if prev is None or rank_key(rec) < rank_key(prev[0]):
                best_by_species[sp] = (rec, ref_name, nm, sub, dmg, adj, as_score)

        if not best_by_species:
            n_unclassified_in_bam += 1
            dec_fh.write(f"{qname}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t"
                        f"UNCLASSIFIED\tno_resolvable_species\n")
            return

        species_hit_count = len(best_by_species)
        nonoryza = sorted(
            ((sp, v) for sp, v in best_by_species.items() if sp not in oryza_set),
            key=lambda kv: rank_key(kv[1][0]))
        oryza = sorted(
            ((sp, v) for sp, v in best_by_species.items() if sp in oryza_set),
            key=lambda kv: rank_key(kv[1][0]))

        rows = []
        for rank, (sp, v) in enumerate(nonoryza[:args.top_n], start=1):
            rows.append((sp, v, rank, False))
        for sp, v in oryza:
            rows.append((sp, v, "NA", True))

        for sp, (rec, ref_name, nm, sub, dmg, adj, as_score), rank, always in rows:
            top10_fh.write(f"{qname}\t{qlen}\t{species_hit_count}\t{rank}\t"
                          f"{'yes' if always else 'no'}\t{sp}\t{tax.sci_name(sp)}\t"
                          f"{ref_name}\t{nm}\t{'NA' if sub is None else sub}\t{dmg}\t"
                          f"{adj}\t{as_score}\n")

        best_nonoryza = nonoryza[0] if nonoryza else None
        best_oryza = oryza[0] if oryza else None

        def cols(entry):
            sp, (rec, ref_name, nm, sub, dmg, adj, as_score) = entry
            return sp, tax.sci_name(sp), nm, dmg, adj

        if best_oryza is None:
            decision, reason = "REJECT", "no_oryza_hit"
            n_rejected_no_oryza += 1
            no_t, no_n, no_nm, no_dmg, no_adj = cols(best_nonoryza)
            o_t = o_n = o_nm = o_dmg = o_adj = "NA"
        elif best_nonoryza is None:
            decision, reason = "KEEP", "oryza_only_no_competitor"
            n_kept += 1
            keep_names.add(qname)
            o_t, o_n, o_nm, o_dmg, o_adj = cols(best_oryza)
            no_t = no_n = no_nm = no_dmg = no_adj = "NA"
        else:
            no_t, no_n, no_nm, no_dmg, no_adj = cols(best_nonoryza)
            o_t, o_n, o_nm, o_dmg, o_adj = cols(best_oryza)
            if o_adj <= no_adj:
                decision, reason = "KEEP", "oryza_at_least_as_good"
                n_kept += 1
                keep_names.add(qname)
            else:
                decision, reason = "REJECT", "nonoryza_better"
                n_rejected_nonoryza_better += 1

        dec_fh.write(f"{qname}\t{no_t}\t{no_n}\t{no_nm}\t{no_dmg}\t{no_adj}\t"
                    f"{o_t}\t{o_n}\t{o_nm}\t{o_dmg}\t{o_adj}\t{decision}\t{reason}\n")

    print("[run] streaming BAM ...", file=sys.stderr)
    cur_qname, cur_qlen, cur_alns = None, 0, []
    n_distinct = 0
    for aln in bam:
        if aln.is_unmapped:
            continue
        qname = aln.query_name
        if cur_qname is not None and qname != cur_qname:
            flush_read(cur_qname, cur_qlen, cur_alns)
            n_distinct += 1
            if args.limit_reads is not None and n_distinct >= args.limit_reads:
                cur_qname = None
                break
            cur_alns = []
            cur_qlen = 0
        if cur_qlen == 0:
            ql = aln.query_length or aln.infer_query_length() or 0
            cur_qlen = ql
        cur_qname = qname
        cur_alns.append(aln)
    if cur_qname is not None:
        flush_read(cur_qname, cur_qlen, cur_alns)
        n_distinct += 1
    bam.close()

    top10_fh.close()
    dec_fh.close()
    os.rename(top10_tmp, top10_tmp[:-4])
    os.rename(dec_tmp, dec_tmp[:-4])

    if n_missing_md_or_seq:
        print(f"[warn] {n_missing_md_or_seq} alignments had no usable SEQ/MD "
              f"(substitution_count/damage could not be computed; adjusted_NM "
              f"fell back to NM for those). Spot-check `samtools view` if this "
              f"is a large fraction.", file=sys.stderr)

    print(f"[run] extracting KEEP reads from candidate FASTQ ...", file=sys.stderr)
    fastq_out = os.path.join(args.outdir, f"{args.sample}.besthit_oryza.fastq.gz")
    n_fastq_out, written = extract_fastq(args.fastq, keep_names, fastq_out, args.threads)
    if n_fastq_out != len(keep_names):
        missing = keep_names - written
        print(f"[warn] {len(missing)} KEEP read names were not found in "
              f"{args.fastq} (id mismatch?). Example: "
              f"{next(iter(missing), None)}", file=sys.stderr)

    if smoke:
        input_reads = n_distinct
    else:
        input_reads = count_fastq_reads(args.fastq)
    unclassified_reads = n_unclassified_in_bam + (input_reads - n_reads_with_alignment)

    check_sum = (n_kept + n_rejected_nonoryza_better + n_rejected_no_oryza +
                n_rejected_low_quality + unclassified_reads)
    consistent = (check_sum == input_reads) and (n_fastq_out == n_kept)

    # Written even on a failed consistency check, so the numbers are there to
    # debug -- only the .finished marker below is gated on success.
    summary_path = os.path.join(args.outdir, f"{args.sample}.summary.tsv")
    summary_tmp = summary_path + ".tmp"
    with open(summary_tmp, "w") as fh:
        fh.write("sample\tinput_reads\treads_with_alignment\treads_with_oryza_hit\t"
                "kept_reads\trejected_nonoryza_better\trejected_no_oryza\t"
                "rejected_low_quality\tunclassified_reads\n")
        fh.write(f"{args.sample}\t{input_reads}\t{n_reads_with_alignment}\t"
                f"{n_kept + n_rejected_nonoryza_better}\t{n_kept}\t"
                f"{n_rejected_nonoryza_better}\t{n_rejected_no_oryza}\t"
                f"{n_rejected_low_quality}\t{unclassified_reads}\n")
    os.rename(summary_tmp, summary_path)

    print(f"[summary] {args.sample}: input={input_reads} "
          f"with_alignment={n_reads_with_alignment} "
          f"oryza_hit={n_kept + n_rejected_nonoryza_better} kept={n_kept} "
          f"rejected_nonoryza_better={n_rejected_nonoryza_better} "
          f"rejected_no_oryza={n_rejected_no_oryza} "
          f"rejected_low_quality={n_rejected_low_quality} "
          f"unclassified={unclassified_reads}", file=sys.stderr)

    if not consistent:
        msg = (f"[check] input_reads={input_reads} != kept({n_kept}) + "
              f"rejected_nonoryza_better({n_rejected_nonoryza_better}) + "
              f"rejected_no_oryza({n_rejected_no_oryza}) + "
              f"rejected_low_quality({n_rejected_low_quality}) + "
              f"unclassified({unclassified_reads}) = {check_sum}; "
              f"fastq_out={n_fastq_out} vs kept={n_kept}")
        if smoke:
            print(f"[smoke-test] {msg} (not enforced in smoke-test mode)",
                  file=sys.stderr)
        else:
            print(f"[error] {msg}", file=sys.stderr)
            sys.exit(1)

    if smoke:
        print("[done] (smoke test -- no .finished marker written)", file=sys.stderr)
    else:
        open(os.path.join(args.outdir, f"{args.sample}.finished"), "w").close()
        print("[done]", file=sys.stderr)


if __name__ == "__main__":
    main()
