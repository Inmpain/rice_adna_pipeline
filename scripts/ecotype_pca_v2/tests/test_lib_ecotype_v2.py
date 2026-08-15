#!/usr/bin/env python3
"""Pure-python unit tests for lib_ecotype_v2.py -- no server, no plink2/
convertf/bedtools, no real panel data. Actually executed (not just written)
as part of the Batch 1 correction pass; see the commit message / chat report
for the pass/fail transcript. Covers item 15's requirements that don't need
external tools: TV/ALL classification, primary/S1-S4 sensitivity resolution
(including invalid-enum rejection), paperlike_5kb window boundary math,
duplicate-ID detection, and Panel A structural validation.

Run with: python3 -m unittest test_lib_ecotype_v2 -v
(from scripts/ecotype_pca_v2/tests/, with scripts/ecotype_pca_v2/ on
PYTHONPATH, e.g. `PYTHONPATH=.. python3 -m unittest test_lib_ecotype_v2 -v`)
"""
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from lib_ecotype_v2 import (is_transversion, resolve_marker_params, genomic_window_index,
                             find_duplicate_ids, validate_panel_a_bim)


def make_config():
    return {
        "version": "ecotype_pca_v2",
        "panel_A_3k": {
            "axis_labels": ["IND", "AUS", "ARO", "TRJ", "TEJ"],
            "project_labels": ["ADM"],
            "geno": 0.05, "maf": 0.01, "ld_window_kb": 100, "ld_r2": 0.20,
            "sensitivity_S4_maf": 0.05,
        },
        "panel_B_720": {
            "axis_mode": "all_modern",
            "geno": 0.10, "maf": 0.01,
            "ld_primary": {"window_kb": 100, "r2": 0.20},
            "paperlike_5kb": {"window_bp": 5000, "seed": 20260814},
        },
        "panel_C_civan": {
            "axis_labels": ["indica", "aus", "aromatic", "japonica",
                             "japonica_temperate", "japonica_tropical"],
            "expected_axis_builder_n": 595,
            "geno": 0.05, "maf": 0.01, "ld_window_kb": 100, "ld_r2": 0.20,
        },
        "sensitivity": [
            {"label": "S1", "window_kb": 50, "r2": 0.20},
            {"label": "S2", "window_kb": 200, "r2": 0.20},
            {"label": "S3", "window_kb": 100, "r2": 0.10},
        ],
    }


class TestTransversion(unittest.TestCase):
    def test_transversions_true(self):
        for a1, a2 in [("A", "C"), ("A", "T"), ("C", "G"), ("G", "T")]:
            self.assertTrue(is_transversion(a1, a2), f"{a1}/{a2} should be TV")
            self.assertTrue(is_transversion(a2, a1), f"{a2}/{a1} should be TV (order-independent)")

    def test_transitions_false(self):
        for a1, a2 in [("A", "G"), ("C", "T")]:
            self.assertFalse(is_transversion(a1, a2), f"{a1}/{a2} should NOT be TV")

    def test_identical_or_malformed_none(self):
        self.assertIsNone(is_transversion("A", "A"))
        self.assertIsNone(is_transversion("AC", "T"))
        self.assertIsNone(is_transversion("A", "N"))

    def test_case_insensitive(self):
        self.assertTrue(is_transversion("a", "c"))


class TestResolveMarkerParams(unittest.TestCase):
    def setUp(self):
        self.cfg = make_config()

    def test_primary_panel_A(self):
        p = resolve_marker_params(self.cfg, "A", "primary")
        self.assertEqual(p["geno"], 0.05)
        self.assertEqual(p["maf"], 0.01)
        self.assertEqual(p["ld_window_kb"], 100)
        self.assertEqual(p["ld_r2"], 0.20)

    def test_S1_S2_S3_shared_grid_applies_to_all_panels(self):
        for panel in ("A", "B", "C"):
            s1 = resolve_marker_params(self.cfg, panel, "S1")
            self.assertEqual((s1["ld_window_kb"], s1["ld_r2"]), (50, 0.20), panel)
            s2 = resolve_marker_params(self.cfg, panel, "S2")
            self.assertEqual((s2["ld_window_kb"], s2["ld_r2"]), (200, 0.20), panel)
            s3 = resolve_marker_params(self.cfg, panel, "S3")
            self.assertEqual((s3["ld_window_kb"], s3["ld_r2"]), (100, 0.10), panel)

    def test_S4_only_valid_for_panel_A(self):
        p = resolve_marker_params(self.cfg, "A", "S4")
        self.assertEqual(p["maf"], 0.05)
        self.assertEqual(p["ld_window_kb"], 100)  # primary LD window unchanged
        for panel in ("B", "C"):
            with self.assertRaises(SystemExit):
                resolve_marker_params(self.cfg, panel, "S4")

    def test_panel_B_uses_ld_primary_block(self):
        p = resolve_marker_params(self.cfg, "B", "primary")
        self.assertEqual(p["ld_window_kb"], 100)
        self.assertEqual(p["ld_r2"], 0.20)
        self.assertEqual(p["geno"], 0.10)

    def test_unknown_panel_letter_hard_fails(self):
        with self.assertRaises(SystemExit):
            resolve_marker_params(self.cfg, "D", "primary")
        with self.assertRaises(SystemExit):
            resolve_marker_params(self.cfg, "a", "primary")  # lowercase typo

    def test_unknown_sensitivity_hard_fails(self):
        with self.assertRaises(SystemExit):
            resolve_marker_params(self.cfg, "A", "S5")
        with self.assertRaises(SystemExit):
            resolve_marker_params(self.cfg, "A", "thinning_only")  # the old broken interface
        with self.assertRaises(SystemExit):
            resolve_marker_params(self.cfg, "A", "Primary")  # case typo


class TestWindowBoundary(unittest.TestCase):
    def test_1based_boundary_exact(self):
        # non-overlapping 5000bp windows over 1-based coordinates:
        # window 0 = positions 1..5000, window 1 = 5001..10000, etc.
        self.assertEqual(genomic_window_index(1, 5000), 0)
        self.assertEqual(genomic_window_index(4999, 5000), 0)
        self.assertEqual(genomic_window_index(5000, 5000), 0)  # the bug this fixes
        self.assertEqual(genomic_window_index(5001, 5000), 1)
        self.assertEqual(genomic_window_index(10000, 5000), 1)
        self.assertEqual(genomic_window_index(10001, 5000), 2)

    def test_other_window_sizes(self):
        self.assertEqual(genomic_window_index(100, 100), 0)
        self.assertEqual(genomic_window_index(101, 100), 1)
        self.assertEqual(genomic_window_index(1, 1), 0)
        self.assertEqual(genomic_window_index(2, 1), 1)


class TestDuplicateIds(unittest.TestCase):
    def test_no_duplicates(self):
        self.assertEqual(find_duplicate_ids(["a", "b", "c"]), {})

    def test_duplicates_found(self):
        dups = find_duplicate_ids(["a", "b", "a", "c", "b", "b"])
        self.assertEqual(dups, {"a": 2, "b": 3})

    def test_empty(self):
        self.assertEqual(find_duplicate_ids([]), {})


class TestValidatePanelABim(unittest.TestCase):
    def _write_bim(self, rows):
        fh = tempfile.NamedTemporaryFile(mode="w", suffix=".bim", delete=False)
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
        fh.close()
        return fh.name

    def test_clean_bim_no_problems(self):
        path = self._write_bim([
            ("1", "snp1", "0", "100", "A", "C"),
            ("1", "snp2", "0", "200", "G", "T"),
            ("12", "snp3", "0", "300", "A", "G"),
        ])
        try:
            self.assertEqual(validate_panel_a_bim(path), [])
        finally:
            os.unlink(path)

    def test_bad_chrom(self):
        path = self._write_bim([
            ("1", "snp1", "0", "100", "A", "C"),
            ("Un", "snp2", "0", "200", "G", "T"),
        ])
        try:
            problems = validate_panel_a_bim(path)
            self.assertTrue(any("chromosomes outside 1-12" in p for p in problems))
        finally:
            os.unlink(path)

    def test_non_biallelic(self):
        path = self._write_bim([
            ("1", "snp1", "0", "100", "A", "A"),   # identical alleles
            ("1", "snp2", "0", "200", "AC", "G"),  # not single-base
        ])
        try:
            problems = validate_panel_a_bim(path)
            self.assertTrue(any("non-biallelic" in p for p in problems))
        finally:
            os.unlink(path)

    def test_duplicate_snp_id(self):
        path = self._write_bim([
            ("1", "snpX", "0", "100", "A", "C"),
            ("1", "snpX", "0", "200", "G", "T"),
        ])
        try:
            problems = validate_panel_a_bim(path)
            self.assertTrue(any("duplicate SNP IDs" in p for p in problems))
        finally:
            os.unlink(path)

    def test_duplicate_position(self):
        path = self._write_bim([
            ("1", "snp1", "0", "100", "A", "C"),
            ("1", "snp2", "0", "100", "G", "T"),
        ])
        try:
            problems = validate_panel_a_bim(path)
            self.assertTrue(any("duplicate (chrom,pos)" in p for p in problems))
        finally:
            os.unlink(path)


if __name__ == "__main__":
    unittest.main()
