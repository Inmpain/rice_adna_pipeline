#!/usr/bin/env python3
"""
Per-read Oryza vs non-Oryza competitive best-hit filter, with ancient-DNA
terminal-damage-aware NM correction.

This does NOT reuse besthit_competitive.py's hierarchical ngsLCA-style walk --
the decision here is a flat, two-sided competition (Oryza vs everything else),
and NM has to be corrected for expected 5'C->T / 3'G->A deamination before the
competition is judged, which the generic classifier does not do.

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
(--top-n); the Oryza species (taxids in --oryza-taxids) are ALWAYS kept in
addition, never competing for one of the 10 non-Oryza slots.

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
    input_reads == kept_reads + rejected_nonoryza_better + rejected_no_oryza + unclassified_reads
holds exactly (checked at the end, hard failure unless --limit-reads is set).

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
                                           never race on one shared file)
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
        self._species_cache = {}
        with open(nodes_file) as fh:
            for line in fh:
                p = line.split("|")
                if len(p) < 3:
                    continue
                tid = p[0].strip()
                self.parent[tid] = p[1].strip()
                self.rank[tid] = p[2].strip()
        with open(names_file) as fh:
            for line in fh:
                if "scientific name" not in line:
                    continue
                p = line.split("|")
                if len(p) < 2:
                    continue
                self.name[p[0].strip()] = p[1].strip()

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
# Per-alignment metrics
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
    ap.add_argument("--oryza-taxids", nargs="+", default=["4529", "4530", "4536"],
                    help="default: O. rufipogon, O. sativa, O. nivara")
    ap.add_argument("--damage-window", type=int, default=5)
    ap.add_argument("--top-n", type=int, default=10,
                    help="non-Oryza species kept in the audit table per read")
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
    for t in args.oryza_taxids:
        if t not in tax.parent:
            sys.exit(f"[error] --oryza-taxids {t} not found in nodes.dmp")
    oryza_set = set(args.oryza_taxids)

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
    n_missing_md_or_seq = 0

    def flush_read(qname, qlen, alns):
        """alns: list of pysam.AlignedSegment for this one read."""
        nonlocal n_reads_with_alignment, n_unclassified_in_bam, n_kept
        nonlocal n_rejected_nonoryza_better, n_rejected_no_oryza, n_missing_md_or_seq

        n_reads_with_alignment += 1
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

    check_sum = n_kept + n_rejected_nonoryza_better + n_rejected_no_oryza + unclassified_reads
    consistent = (check_sum == input_reads) and (n_fastq_out == n_kept)

    # Written even on a failed consistency check, so the numbers are there to
    # debug -- only the .finished marker below is gated on success.
    summary_path = os.path.join(args.outdir, f"{args.sample}.summary.tsv")
    summary_tmp = summary_path + ".tmp"
    with open(summary_tmp, "w") as fh:
        fh.write("sample\tinput_reads\treads_with_alignment\treads_with_oryza_hit\t"
                "kept_reads\trejected_nonoryza_better\trejected_no_oryza\t"
                "unclassified_reads\n")
        fh.write(f"{args.sample}\t{input_reads}\t{n_reads_with_alignment}\t"
                f"{n_kept + n_rejected_nonoryza_better}\t{n_kept}\t"
                f"{n_rejected_nonoryza_better}\t{n_rejected_no_oryza}\t"
                f"{unclassified_reads}\n")
    os.rename(summary_tmp, summary_path)

    print(f"[summary] {args.sample}: input={input_reads} "
          f"with_alignment={n_reads_with_alignment} "
          f"oryza_hit={n_kept + n_rejected_nonoryza_better} kept={n_kept} "
          f"rejected_nonoryza_better={n_rejected_nonoryza_better} "
          f"rejected_no_oryza={n_rejected_no_oryza} "
          f"unclassified={unclassified_reads}", file=sys.stderr)

    if not consistent:
        msg = (f"[check] input_reads={input_reads} != kept({n_kept}) + "
              f"rejected_nonoryza_better({n_rejected_nonoryza_better}) + "
              f"rejected_no_oryza({n_rejected_no_oryza}) + "
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
