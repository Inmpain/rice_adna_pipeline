import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).parents[3]


class Contract(unittest.TestCase):
    def test_scripts_help(self):
        names = [
            "09_export_fixed_reference_eigenstrat.py",
            "10_call_ancient_fixed_markers.py",
            "11_build_ancient_callability.py",
            "12_build_ancient_overlap_matrix.py",
            "13_merge_ancients_fixed_panel.py",
            "15_pca_qc.py",
            "16_projection_summary.py",
            "17_exact_mask_validation.py",
            "18_validation_metrics.py",
        ]
        for name in names:
            result = subprocess.run(
                [sys.executable, str(ROOT / "scripts/ecotype_pca_v2" / name), "--help"],
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, name)
        result = subprocess.run(
            ["bash", str(ROOT / "scripts/ecotype_pca_v2/14_run_fixed_smartpca.sh"), "--help"]
        )
        self.assertEqual(result.returncode, 0)

    def test_shared_seed_and_encoding(self):
        sys.path.insert(0, str(ROOT / "scripts/ecotype_pca_v2"))
        from fixed_projection_lib import VALID_GENOTYPES, stable_int_seed

        self.assertEqual(stable_int_seed("x", 1), stable_int_seed("x", 1))
        self.assertEqual(VALID_GENOTYPES, frozenset("0129"))

    def test_pca_qc_allows_only_missing_projection_samples(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            ind = root / "merged.ind"
            evec = root / "pca.evec"
            out = root / "qc.tsv"
            ind.write_text("modern1\tU\tindica\nancient0\tU\tAncient\n")
            evec.write_text("#eigvals\nmodern1 " + " ".join(["0.1"] * 10) + " indica\n")
            command = [
                sys.executable,
                str(ROOT / "scripts/ecotype_pca_v2/15_pca_qc.py"),
                "--evec",
                str(evec),
                "--ind",
                str(ind),
                "--expected-n",
                "2",
                "--out",
                str(out),
            ]

            ok = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(ok.returncode, 0, ok.stderr)
            with out.open(newline="") as handle:
                qc = {
                    row["metric"]: row["value"]
                    for row in csv.DictReader(handle, delimiter="\t")
                }
            self.assertEqual(qc["missing_sample_ids"], "ancient0")

            out.unlink()
            ind.write_text("modern0\tU\tindica\nancient1\tU\tAncient\n")
            evec.write_text("#eigvals\nancient1 " + " ".join(["0.1"] * 10) + " Ancient\n")
            bad = subprocess.run(command, capture_output=True, text=True)
            self.assertNotEqual(bad.returncode, 0)
            self.assertIn("omitted non-projection", bad.stderr)


if __name__ == "__main__":
    unittest.main()
