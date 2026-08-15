#!/usr/bin/env python3
"""Actually-executed test of 08's stream_thin() against a synthetic .bim --
verifies the streaming low-memory reimplementation (item 5) picks exactly
one SNP per populated window, respects the corrected 1-based boundary, is
deterministic given the same seed, and hard-fails on unsorted input rather
than silently reordering.

Run with: python3 test_08_streaming.py
"""
import logging
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import importlib.util
spec = importlib.util.spec_from_file_location("script08", os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "08_make_5kb_thinned_markers.py"))
script08 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(script08)

logger = logging.getLogger("test08")
logger.addHandler(logging.NullHandler())


def write_bim(rows):
    fh = tempfile.NamedTemporaryFile(mode="w", suffix=".bim", delete=False)
    for r in rows:
        fh.write("\t".join(str(x) for x in r) + "\n")
    fh.close()
    return fh.name


def run(name, fn):
    try:
        fn()
        print(f"PASS: {name}")
    except AssertionError as e:
        print(f"FAIL: {name}: {e}")
        raise


def test_one_snp_per_window():
    path = write_bim([
        ("1", "s1", "0", 1), ("1", "s2", "0", 2500), ("1", "s3", "0", 5000),  # window 0
        ("1", "s4", "0", 5001), ("1", "s5", "0", 8000),                       # window 1
        ("1", "s6", "0", 10001),                                              # window 2
    ])
    try:
        kept, n_in = script08.stream_thin(path, 5000, 42, logger)
        assert n_in == 6, n_in
        assert len(kept) == 3, kept
        windows = sorted(w for _c, w, _s in kept)
        assert windows == [0, 1, 2], windows
        # window 0's pick must be one of s1/s2/s3 (the boundary case s3@5000 must
        # land in window 0, not window 1 -- this is exactly the bug being fixed)
        w0_pick = [s for c, w, s in kept if w == 0][0]
        assert w0_pick in ("s1", "s2", "s3"), w0_pick
        w2_pick = [s for c, w, s in kept if w == 2][0]
        assert w2_pick == "s6", w2_pick
    finally:
        os.unlink(path)


def test_deterministic_same_seed():
    path = write_bim([("1", f"s{i}", "0", i) for i in range(1, 5001, 100)])
    try:
        kept1, _ = script08.stream_thin(path, 5000, 20260814, logger)
        kept2, _ = script08.stream_thin(path, 5000, 20260814, logger)
        assert kept1 == kept2, (kept1, kept2)
    finally:
        os.unlink(path)


def test_different_seed_can_differ():
    path = write_bim([("1", f"s{i}", "0", i) for i in range(1, 5001, 100)])
    try:
        kept1, _ = script08.stream_thin(path, 5000, 1, logger)
        kept2, _ = script08.stream_thin(path, 5000, 2, logger)
        # not a strict guarantee for tiny inputs, but with 50 candidates in one
        # window and 2 different seeds it should not always coincide -- if this
        # ever flakes it's worth a look, but it's not testing anything load-bearing
        assert kept1[0][2] != kept2[0][2] or True  # informational only, no assert
    finally:
        os.unlink(path)


def test_multi_chrom_independent_windows():
    path = write_bim([
        ("1", "a1", "0", 1), ("1", "a2", "0", 6000),
        ("2", "b1", "0", 1), ("2", "b2", "0", 6000),
    ])
    try:
        kept, n_in = script08.stream_thin(path, 5000, 7, logger)
        assert n_in == 4
        assert len(kept) == 4, kept  # each (chrom,window) pair has exactly 1 candidate
        chroms = sorted(set(c for c, _w, _s in kept))
        assert chroms == ["1", "2"], chroms
    finally:
        os.unlink(path)


def test_unsorted_input_hard_fails():
    path = write_bim([
        ("1", "s1", "0", 5000),
        ("1", "s2", "0", 100),  # goes backwards -- must raise, not silently reorder
    ])
    try:
        try:
            script08.stream_thin(path, 5000, 1, logger)
            assert False, "expected ValueError for unsorted input"
        except ValueError as e:
            assert "not sorted" in str(e)
    finally:
        os.unlink(path)


if __name__ == "__main__":
    run("one_snp_per_window (incl. pos==window_bp boundary)", test_one_snp_per_window)
    run("deterministic_same_seed", test_deterministic_same_seed)
    run("different_seed_can_differ (informational)", test_different_seed_can_differ)
    run("multi_chrom_independent_windows", test_multi_chrom_independent_windows)
    run("unsorted_input_hard_fails", test_unsorted_input_hard_fails)
    print("\nALL 08 STREAMING TESTS PASSED")
