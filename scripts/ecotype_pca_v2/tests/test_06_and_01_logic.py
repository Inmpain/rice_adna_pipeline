#!/usr/bin/env python3
"""Actually-executed tests for 06's per-panel reference-set builders (missing
axis label / count-mismatch hard fail, item 11) and 01's duplicate-ID
detection (item 12), using tiny synthetic fixtures -- no plink2/bedtools/
convertf, no server. Complements test_lib_ecotype_v2.py and
test_08_streaming.py for item 15's "duplicate ID / missing label / wrong
enum" coverage requirement.

Run with: python3 test_06_and_01_logic.py
"""
import importlib.util
import logging
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
PARENT = os.path.join(HERE, "..")


def load_module(name, filename):
    spec = importlib.util.spec_from_file_location(name, os.path.join(PARENT, filename))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


script06 = load_module("script06", "06_build_reference_sample_set.py")
script01 = load_module("script01", "01_make_panel_manifest.py")

logger = logging.getLogger("test")
logger.addHandler(logging.NullHandler())


def make_config():
    return {
        "version": "ecotype_pca_v2",
        "panel_A_3k": {"axis_labels": ["IND", "AUS", "ARO", "TRJ", "TEJ"],
                        "project_labels": ["ADM"]},
        "panel_C_civan": {"axis_labels": ["indica", "aus", "aromatic", "japonica",
                                            "japonica_temperate", "japonica_tropical"],
                            "expected_axis_builder_n": 6},
    }


def rows(*label_counts):
    """label_counts: [(label, n), ...] -> list of {'id','sex','label'} dicts"""
    out = []
    i = 0
    for label, n in label_counts:
        for _ in range(n):
            i += 1
            out.append({"id": f"s{i}", "sex": "U", "label": label})
    return out


def run(name, fn):
    fn()
    print(f"PASS: {name}")


def test_panel_A_missing_label_hard_fails():
    cfg = make_config()
    r = rows(("IND", 3), ("AUS", 2), ("ARO", 1), ("TRJ", 1))  # TEJ missing entirely
    try:
        script06.build_panel_A(cfg, r, logger)
        assert False, "expected SystemExit for missing TEJ label"
    except SystemExit:
        pass


def test_panel_A_all_five_present_ok():
    cfg = make_config()
    r = rows(("IND", 3), ("AUS", 2), ("ARO", 1), ("TRJ", 1), ("TEJ", 4), ("ADM", 2), ("Ancient", 1))
    keep, other = script06.build_panel_A(cfg, r, logger)
    assert len(keep) == 3 + 2 + 1 + 1 + 4, len(keep)
    assert other["project:ADM"] == 2
    assert other["unclassified:Ancient"] == 1


def test_panel_C_count_mismatch_hard_fails():
    cfg = make_config()  # expected_axis_builder_n=6 for this synthetic config
    r = rows(("indica", 1), ("aus", 1), ("aromatic", 1), ("japonica", 1),
              ("japonica_temperate", 1), ("japonica_tropical", 2))  # total 7, not 6
    try:
        script06.build_panel_C(cfg, r, logger)
        assert False, "expected SystemExit for count mismatch (7 != 6)"
    except SystemExit:
        pass


def test_panel_C_exact_count_ok():
    cfg = make_config()
    r = rows(("indica", 1), ("aus", 1), ("aromatic", 1), ("japonica", 1),
              ("japonica_temperate", 1), ("japonica_tropical", 1))  # total 6 == expected
    keep, other = script06.build_panel_C(cfg, r, logger)
    assert len(keep) == 6


def test_panel_C_label_spelling_mismatch_hard_fails():
    """Directly exercises the real Civan-label-parentheses risk flagged in the
    Batch 1 report: if the .ind file actually has 'japonica_(temperate)' but
    config says 'japonica_temperate', that label has zero matches -- must
    hard-fail, never silently proceed with 5 of 6 groups."""
    cfg = make_config()
    r = rows(("indica", 1), ("aus", 1), ("aromatic", 1), ("japonica", 1),
              ("japonica_(temperate)", 1), ("japonica_(tropical)", 1))  # parens, doesn't match config
    try:
        script06.build_panel_C(cfg, r, logger)
        assert False, "expected SystemExit for label spelling mismatch"
    except SystemExit:
        pass


def write_ind(rows_):
    fh = tempfile.NamedTemporaryFile(mode="w", suffix=".ind", delete=False)
    for r in rows_:
        fh.write(f"{r['id']}\t{r['sex']}\t{r['label']}\n")
    fh.close()
    return fh.name


def write_snp(rows_):
    fh = tempfile.NamedTemporaryFile(mode="w", suffix=".snp", delete=False)
    for snpid, chrom, pos in rows_:
        fh.write(f"{snpid}\t{chrom}\t0\t{pos}\tA\tC\n")
    fh.close()
    return fh.name


def test_01_duplicate_sample_id_fails():
    ind_path = write_ind([{"id": "s1", "sex": "U", "label": "IND"},
                           {"id": "s1", "sex": "U", "label": "AUS"}])
    snp_path = write_snp([("snp1", "1", 100), ("snp2", "1", 200)])
    d = tempfile.mkdtemp()
    prefix = "panel"
    os.rename(ind_path, os.path.join(d, prefix + ".ind"))
    os.rename(snp_path, os.path.join(d, prefix + ".snp"))
    cfg = {"inputs": {"panel_A_3k": {"dir": d, "prefix": prefix, "filtered_suffix": ".filtered"}}}
    row, ok = script01.manifest_for_panel(cfg, "panel_A_3k", logger)
    assert ok is False, "expected ok=False for duplicate sample ID"


def test_01_duplicate_snp_id_fails():
    d = tempfile.mkdtemp()
    prefix = "panel"
    with open(os.path.join(d, prefix + ".ind"), "w") as fh:
        fh.write("s1\tU\tIND\n")
    with open(os.path.join(d, prefix + ".snp"), "w") as fh:
        fh.write("snpX\t1\t0\t100\tA\tC\n")
        fh.write("snpX\t1\t0\t200\tG\tT\n")
    cfg = {"inputs": {"panel_A_3k": {"dir": d, "prefix": prefix, "filtered_suffix": ".filtered"}}}
    row, ok = script01.manifest_for_panel(cfg, "panel_A_3k", logger)
    assert ok is False, "expected ok=False for duplicate SNP ID"


def test_01_clean_input_ok():
    d = tempfile.mkdtemp()
    prefix = "panel"
    with open(os.path.join(d, prefix + ".ind"), "w") as fh:
        fh.write("s1\tU\tIND\ns2\tU\tAUS\n")
    with open(os.path.join(d, prefix + ".snp"), "w") as fh:
        fh.write("snp1\t1\t0\t100\tA\tC\nsnp2\t1\t0\t200\tG\tT\n")
    cfg = {"inputs": {"panel_A_3k": {"dir": d, "prefix": prefix, "filtered_suffix": ".filtered"}}}
    row, ok = script01.manifest_for_panel(cfg, "panel_A_3k", logger)
    assert ok is True, "expected ok=True for clean input"
    assert row["raw_n_samples"] == 2
    assert row["n_snps"] == 2


if __name__ == "__main__":
    run("panel_A_missing_label_hard_fails", test_panel_A_missing_label_hard_fails)
    run("panel_A_all_five_present_ok", test_panel_A_all_five_present_ok)
    run("panel_C_count_mismatch_hard_fails", test_panel_C_count_mismatch_hard_fails)
    run("panel_C_exact_count_ok", test_panel_C_exact_count_ok)
    run("panel_C_label_spelling_mismatch_hard_fails (the real Civan parens risk)",
        test_panel_C_label_spelling_mismatch_hard_fails)
    run("01_duplicate_sample_id_fails", test_01_duplicate_sample_id_fails)
    run("01_duplicate_snp_id_fails", test_01_duplicate_snp_id_fails)
    run("01_clean_input_ok", test_01_clean_input_ok)
    print("\nALL 06/01 LOGIC TESTS PASSED")
