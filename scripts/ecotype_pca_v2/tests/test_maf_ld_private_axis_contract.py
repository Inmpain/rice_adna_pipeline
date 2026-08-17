import subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).parents[3]


class Contract(unittest.TestCase):
    def test_scripts_help(self):
        for name in ['23_validate_snp_ref_against_fasta.py',
                     '24_extract_sample_covered_sites.py',
                     '25_intersect_snplists.py',
                     '27_ancient_coverage_first_ld_prune.py']:
            r = subprocess.run(
                [sys.executable, str(ROOT / 'scripts/ecotype_pca_v2' / name), '--help'],
                capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, f"{name}: {r.stderr}")

    def test_bash_help_and_syntax(self):
        for name in ['26_plot_pc_pairs.sh']:
            r = subprocess.run(
                ['bash', str(ROOT / 'scripts/ecotype_pca_v2' / name), '--help'],
                capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, f"{name}: {r.stderr}")
        for name in ['07_make_fixed_markers.sh',
                      'workflow/runners/51_civan_maf_ld_and_private_axis.sh']:
            r = subprocess.run(
                ['bash', '-n', str(ROOT / 'scripts/ecotype_pca_v2' / name)],
                capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, f"{name}: {r.stderr}")

    def test_07_help(self):
        r = subprocess.run(
            ['bash', str(ROOT / 'scripts/ecotype_pca_v2/07_make_fixed_markers.sh'), '--help'],
            capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn('pooled_mixed', r.stdout)


if __name__ == '__main__':
    unittest.main()
