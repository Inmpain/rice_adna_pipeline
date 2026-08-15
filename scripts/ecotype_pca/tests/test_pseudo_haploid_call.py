#!/usr/bin/env python3
"""Actually-executed regression tests for pseudo_haploid_call.py, built
around a real synthetic BAM (via pysam, not a mock of the pileup logic).
No server, no real ancient DNA data.

Core invariant this protects: TV and ALL must make the IDENTICAL
pseudo-haploid call at any site both tracks call. The pre-2026-08-15
version violated this (GPT review of that version, reproduced during the
fix: same seed, same synthetic reads at a shared transversion site,
TV called genotype 0 and ALL called genotype 2 -- purely because ALL had
consumed one extra random.choice() draw at an earlier transition site).
Fixed by deriving each site's RNG only from (seed, contig, position). See
the script's own docstring point 5 and the fix commit message for the
before/after transcript.

Run with: python3 test_pseudo_haploid_call.py
Requires pysam (pip install pysam) and pseudo_haploid_call.py alongside
this file.
"""
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import pysam

FIXED = os.path.join(HERE, "..", "pseudo_haploid_call.py")

PASS = 0
FAIL = 0


def check(name, cond, detail=""):
    global PASS, FAIL
    if cond:
        print(f"PASS: {name}")
        PASS += 1
    else:
        print(f"FAIL: {name}  {detail}")
        FAIL += 1


def build_synthetic_bam(workdir, reads_at_10, reads_at_20, contig="chr01"):
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": contig, "LN": 50}]}
    bam_path = os.path.join(workdir, "synthetic.bam")
    sorted_path = os.path.join(workdir, "synthetic.sorted.bam")
    with pysam.AlignmentFile(bam_path, "wb", header=header) as bam:
        i = 0
        for base, refpos0 in [(b, 9) for b in reads_at_10] + [(b, 19) for b in reads_at_20]:
            i += 1
            a = pysam.AlignedSegment()
            a.query_name = f"read_{i}"
            a.query_sequence = base
            a.flag = 0
            a.reference_id = 0
            a.reference_start = refpos0
            a.mapping_quality = 40
            a.cigar = [(0, 1)]
            a.query_qualities = pysam.qualitystring_to_array("I")
            bam.write(a)
    pysam.sort("-o", sorted_path, bam_path)
    pysam.index(sorted_path)
    return sorted_path


def write_snp(workdir, rows, name="panel.snp"):
    path = os.path.join(workdir, name)
    with open(path, "w") as fh:
        for row in rows:
            fh.write("\t".join(str(x) for x in row) + "\n")
    return path


def run_script(script, args):
    proc = subprocess.run([sys.executable, script] + args, capture_output=True, text=True)
    return proc


def test_fixed_tv_all_agree():
    """The actual fix: same scenario, fixed script -- TV and ALL must now
    agree at every site both tracks call."""
    with tempfile.TemporaryDirectory() as d:
        bam = build_synthetic_bam(d, ["A", "G", "A", "G"], ["A", "C", "A", "C", "A"])
        snp = write_snp(d, [
            ("snp_transition", 1, 0, 10, "A", "G"),
            ("snp_transversion", 1, 0, 20, "A", "C"),
        ])
        tv_out = os.path.join(d, "tv.out")
        all_out = os.path.join(d, "all.out")
        p1 = run_script(FIXED, ["--bam", bam, "--panel-snp", snp, "--out", tv_out, "--seed", "1"])
        p2 = run_script(FIXED, ["--bam", bam, "--panel-snp", snp, "--out", all_out,
                                 "--seed", "1", "--no-transversions-only"])
        check("fixed_script_runs_ok", p1.returncode == 0 and p2.returncode == 0,
              f"tv_rc={p1.returncode} all_rc={p2.returncode}\n{p1.stderr}\n{p2.stderr}")
        tv_lines = open(tv_out).read().split()
        all_lines = open(all_out).read().split()
        check("fixed_TV_and_ALL_agree_at_shared_transversion_site",
              tv_lines[1] == all_lines[1],
              f"expected agreement, got TV={tv_lines[1]} ALL={all_lines[1]}")
        check("fixed_TV_still_skips_transition_site", tv_lines[0] == "9", tv_lines[0])
        check("fixed_ALL_calls_transition_site", all_lines[0] != "9", all_lines[0])


def test_fixed_deterministic_across_repeated_runs():
    """Same seed, same BAM, same panel, run 3 times -> identical output each time."""
    with tempfile.TemporaryDirectory() as d:
        bam = build_synthetic_bam(d, ["A", "G", "T", "G", "C"], ["A", "C", "T", "C", "A", "C"])
        snp = write_snp(d, [
            ("s1", 1, 0, 10, "A", "G"),
            ("s2", 1, 0, 20, "A", "C"),
        ])
        outs = []
        for i in range(3):
            out = os.path.join(d, f"run{i}.out")
            p = run_script(FIXED, ["--bam", bam, "--panel-snp", snp, "--out", out,
                                    "--seed", "42", "--no-transversions-only"])
            check(f"deterministic_run{i}_ok", p.returncode == 0, p.stderr)
            outs.append(open(out).read())
        check("deterministic_all_three_runs_identical", outs[0] == outs[1] == outs[2], outs)


def test_contig_format_mismatch_hard_fails():
    with tempfile.TemporaryDirectory() as d:
        bam = build_synthetic_bam(d, ["A"], ["A"], contig="1")  # named "1", not "chr01"
        snp = write_snp(d, [("s1", 1, 0, 10, "A", "C")])
        out = os.path.join(d, "out.txt")
        p = run_script(FIXED, ["--bam", bam, "--panel-snp", snp, "--out", out])
        check("contig_mismatch_hard_fails_nonzero", p.returncode != 0, p.stderr)
        check("contig_mismatch_error_message_clear", "FATAL" in p.stderr and "contig" in p.stderr.lower())


def test_duplicate_snp_id_hard_fails():
    with tempfile.TemporaryDirectory() as d:
        bam = build_synthetic_bam(d, ["A"], ["A"])
        snp = write_snp(d, [
            ("dup", 1, 0, 10, "A", "C"),
            ("dup", 1, 0, 20, "A", "G"),
        ])
        out = os.path.join(d, "out.txt")
        p = run_script(FIXED, ["--bam", bam, "--panel-snp", snp, "--out", out])
        check("duplicate_snp_id_hard_fails", p.returncode != 0, p.stderr)


def test_invalid_ref_alt_hard_fails():
    with tempfile.TemporaryDirectory() as d:
        bam = build_synthetic_bam(d, ["A"], ["A"])
        snp = write_snp(d, [("s1", 1, 0, 10, "A", "A")])  # ref==alt, invalid
        out = os.path.join(d, "out.txt")
        p = run_script(FIXED, ["--bam", bam, "--panel-snp", snp, "--out", out])
        check("invalid_ref_alt_hard_fails", p.returncode != 0, p.stderr)


def test_missing_bai_hard_fails():
    with tempfile.TemporaryDirectory() as d:
        bam = build_synthetic_bam(d, ["A"], ["A"])
        os.unlink(bam + ".bai")
        snp = write_snp(d, [("s1", 1, 0, 10, "A", "C")])
        out = os.path.join(d, "out.txt")
        p = run_script(FIXED, ["--bam", bam, "--panel-snp", snp, "--out", out])
        check("missing_bai_hard_fails", p.returncode != 0, p.stderr)


def test_malformed_snp_line_hard_fails_not_crashes():
    """A too-short line must produce a clear FATAL, not a raw IndexError
    traceback (the field-length-check-too-late bug)."""
    with tempfile.TemporaryDirectory() as d:
        bam = build_synthetic_bam(d, ["A"], ["A"])
        snp = write_snp(d, [("s1", 1)])  # only 2 fields
        out = os.path.join(d, "out.txt")
        p = run_script(FIXED, ["--bam", bam, "--panel-snp", snp, "--out", out])
        check("malformed_line_hard_fails", p.returncode != 0, p.stderr)
        check("malformed_line_no_raw_traceback", "Traceback" not in p.stderr, p.stderr)


def test_unsupported_chromosome_counted():
    with tempfile.TemporaryDirectory() as d:
        bam = build_synthetic_bam(d, ["A", "C"], ["A", "C"])
        snp = write_snp(d, [
            ("s1", "X", 0, 10, "A", "C"),  # non-numeric chromosome
            ("s2", 1, 0, 20, "A", "C"),
        ])
        out = os.path.join(d, "out.txt")
        report = os.path.join(d, "report.tsv")
        p = run_script(FIXED, ["--bam", bam, "--panel-snp", snp, "--out", out, "--report", report])
        check("unsupported_chrom_run_ok", p.returncode == 0, p.stderr)
        report_text = open(report).read()
        check("unsupported_chromosome_field_present", "unsupported_chromosome\t1" in report_text,
              report_text)
        check("eligible_site_call_rate_field_present",
              "eligible_site_call_rate\t" in report_text, report_text)
        check("allele_match_rate_among_covered_field_present",
              "allele_match_rate_among_covered\t" in report_text, report_text)


def test_two_call_rate_metrics_diverge_as_expected():
    """2026-08-15 correction (GPT review of a4fb1e6): a single 'call_rate'
    conflated two different questions. Construct a scenario with mostly
    no_coverage sites (drags eligible_site_call_rate down) but where every
    site that DID get a read drawn matched a known allele (so
    allele_match_rate_among_covered should be 1.0 regardless) -- the two
    metrics must diverge, proving they're not the same number under a
    different name."""
    with tempfile.TemporaryDirectory() as d:
        bam = build_synthetic_bam(d, ["A"], ["A"])  # only pos 10 and pos 20 covered
        snp_rows = [("covered1", 1, 0, 10, "A", "C"), ("covered2", 1, 0, 20, "A", "C")]
        # add 8 more transversion SNPs with no coverage at all
        for i, pos in enumerate(range(100, 900, 100), start=1):
            snp_rows.append((f"uncovered{i}", 1, 0, pos, "A", "C"))
        snp = write_snp(d, snp_rows)
        out = os.path.join(d, "out.txt")
        report = os.path.join(d, "report.tsv")
        p = run_script(FIXED, ["--bam", bam, "--panel-snp", snp, "--out", out, "--report", report])
        check("diverge_scenario_run_ok", p.returncode == 0, p.stderr)
        values = dict(line.split("\t") for line in open(report).read().splitlines()[1:])
        eligible = float(values["eligible_site_call_rate"])
        match_rate = float(values["allele_match_rate_among_covered"])
        check("eligible_rate_diluted_by_no_coverage", eligible < 0.3, eligible)  # 2/10
        check("match_rate_not_diluted_by_no_coverage", match_rate == 1.0, match_rate)
        check("the_two_metrics_actually_differ", abs(eligible - match_rate) > 0.5,
              (eligible, match_rate))


if __name__ == "__main__":
    test_fixed_tv_all_agree()
    test_fixed_deterministic_across_repeated_runs()
    test_contig_format_mismatch_hard_fails()
    test_duplicate_snp_id_hard_fails()
    test_invalid_ref_alt_hard_fails()
    test_missing_bai_hard_fails()
    test_malformed_snp_line_hard_fails_not_crashes()
    test_unsupported_chromosome_counted()
    test_two_call_rate_metrics_diverge_as_expected()
    print(f"\nTOTAL: {PASS} passed, {FAIL} failed")
    sys.exit(1 if FAIL else 0)
