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

    def test_one_command_server_entry_shell_syntax_and_help(self):
        script = ROOT / 'scripts/ecotype_pca_v2/workflow/submit_coverage_aware_stage50.sh'
        syntax = subprocess.run(['bash', '-n', str(script)], capture_output=True, text=True)
        self.assertEqual(syntax.returncode, 0, syntax.stderr)
        help_result = subprocess.run(['bash', str(script), '--help'], capture_output=True, text=True)
        self.assertEqual(help_result.returncode, 0, help_result.stderr)
        text = script.read_text()
        self.assertNotIn('export CIVAN_UNION_SITES="/path/to/', text)
        self.assertNotIn('export CIVAN_UNION_SITES_TV="/path/to/', text)
        self.assertNotIn('export CIVAN_REFERENCE_KEEP="/path/to/', text)
        self.assertIn('ecotype_pca_workflow.py', text)

    def test_runner_is_fail_closed_and_classifies_every_requested_sample(self):
        runner = (ROOT / 'scripts/ecotype_pca_v2/workflow/runners/50_civan_coverage_aware_projection.sh').read_text()
        self.assertIn('REPORT_ARGS+=("$SAMPLE=', runner)
        self.assertIn('TECHNICAL_FAILURE_N', runner)
        self.assertIn('refusing an incomplete Stage 50 receipt', runner)

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
