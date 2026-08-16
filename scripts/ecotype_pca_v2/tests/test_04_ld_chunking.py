#!/usr/bin/env python3
"""Pure-Python regression tests for 04_audit_720_ld.py chunk geometry."""

import importlib.util
from pathlib import Path


HERE = Path(__file__).resolve().parent
SCRIPT = HERE.parent / "04_audit_720_ld.py"
SPEC = importlib.util.spec_from_file_location("audit_720_ld", SCRIPT)
MOD = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MOD)


def check(name, condition):
    if not condition:
        raise AssertionError(name)
    print(f"PASS: {name}")


def main():
    # Block 1 owns the lower SNP at 19.9Mb.  Its query must include the upper
    # SNP at 20.1Mb across the 20Mb boundary.
    q_from, q_to = MOD.chunk_query_bounds(1, 20_000_000, 45_000_000, 500_000)
    check("first_block_starts_at_owned_range", q_from == 1)
    check("forward_halo_crosses_block_end", q_to == 20_500_000)
    lower, upper = 19_900_000, 20_100_000
    check("cross_boundary_pair_visible_in_owner_block", q_from <= lower <= upper <= q_to)
    check("lower_coordinate_owned_by_first_block", 1 <= min(lower, upper) <= 20_000_000)

    q_from2, q_to2 = MOD.chunk_query_bounds(
        20_000_001, 40_000_000, 40_200_000, 500_000
    )
    check("second_block_has_no_backward_halo", q_from2 == 20_000_001)
    check("forward_halo_clamped_to_chromosome", q_to2 == 40_200_000)

    check("exact_500kb_pair_is_binned", MOD.bin_label(500_000) == "200-500kb")
    check("beyond_window_is_not_binned", MOD.bin_label(500_001) is None)
    print("\nALL 04 LD CHUNKING TESTS PASSED")


if __name__ == "__main__":
    main()
