import subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).parents[3]


class Contract(unittest.TestCase):
    def test_script_help(self):
        r = subprocess.run(
            [sys.executable, str(ROOT / 'scripts/ecotype_pca_v2' / '19_survey_ancient_coverage.py'), '--help'],
            capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)

    def test_format_contig(self):
        sys.path.insert(0, str(ROOT / 'scripts/ecotype_pca_v2'))
        from fixed_projection_lib import format_contig, tally_coverage
        self.assertEqual(format_contig('chr%02d', '1'), 'chr01')
        self.assertEqual(format_contig('chr%02d', '01'), 'chr01')
        self.assertEqual(format_contig('chr%02d', '12'), 'chr12')
        with self.assertRaises(ValueError):
            format_contig('chr%02d', 'X')

    def test_tally_coverage(self):
        sys.path.insert(0, str(ROOT / 'scripts/ecotype_pca_v2'))
        from fixed_projection_lib import tally_coverage
        tally = tally_coverage({'s1': {('chr01', 5), ('chr01', 9)}, 's2': {('chr01', 5)}})
        self.assertEqual(tally[('chr01', 5)], {'s1', 's2'})
        self.assertEqual(tally[('chr01', 9)], {'s1'})


if __name__ == '__main__':
    unittest.main()
