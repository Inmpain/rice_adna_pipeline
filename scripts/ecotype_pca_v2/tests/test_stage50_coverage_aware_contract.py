import subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).parents[3]


class Contract(unittest.TestCase):
    def test_scripts_help(self):
        for name in ['21_extract_fixed_snplist.py', '22_classify_scientific_projection.py']:
            r = subprocess.run(
                [sys.executable, str(ROOT / 'scripts/ecotype_pca_v2' / name), '--help'],
                capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, f"{name}: {r.stderr}")

    def test_runner_shell_syntax(self):
        r = subprocess.run(
            ['bash', '-n', str(ROOT / 'scripts/ecotype_pca_v2/workflow/runners/50_civan_coverage_aware_projection.sh')],
            capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)

    def test_classify_tiers(self):
        sys.path.insert(0, str(ROOT / 'scripts/ecotype_pca_v2'))
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            'classify_mod', ROOT / 'scripts/ecotype_pca_v2' / '22_classify_scientific_projection.py')
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        self.assertEqual(mod.classify(200, 200, 50), 'formal_validation_candidate')
        self.assertEqual(mod.classify(199, 200, 50), 'exploratory_projection')
        self.assertEqual(mod.classify(50, 200, 50), 'exploratory_projection')
        self.assertEqual(mod.classify(49, 200, 50), 'descriptive_only')
        self.assertEqual(mod.classify(0, 200, 50), 'descriptive_only')


if __name__ == '__main__':
    unittest.main()
