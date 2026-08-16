import subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).parents[3]


class Contract(unittest.TestCase):
    def test_script_help(self):
        r = subprocess.run(
            [sys.executable, str(ROOT / 'scripts/ecotype_pca_v2' / '20_filter_coverage_sites_to_transversions.py'), '--help'],
            capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)

    def test_survey_script_help_still_works_after_shared_helper_move(self):
        r = subprocess.run(
            [sys.executable, str(ROOT / 'scripts/ecotype_pca_v2' / '19_survey_ancient_coverage.py'), '--help'],
            capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)


if __name__ == '__main__':
    unittest.main()
