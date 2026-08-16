import subprocess,sys,unittest
from pathlib import Path
ROOT=Path(__file__).parents[3]
class Contract(unittest.TestCase):
 def test_scripts_help(self):
  names=['09_export_fixed_reference_eigenstrat.py','10_call_ancient_fixed_markers.py','11_build_ancient_callability.py','12_build_ancient_overlap_matrix.py','13_merge_ancients_fixed_panel.py','15_pca_qc.py','16_projection_summary.py','17_exact_mask_validation.py','18_validation_metrics.py']
  for n in names:
   r=subprocess.run([sys.executable,str(ROOT/'scripts/ecotype_pca_v2'/n),'--help'],capture_output=True,text=True); self.assertEqual(r.returncode,0,n)
  r=subprocess.run(['bash',str(ROOT/'scripts/ecotype_pca_v2/14_run_fixed_smartpca.sh'),'--help']); self.assertEqual(r.returncode,0)
 def test_shared_seed_and_encoding(self):
  sys.path.insert(0,str(ROOT/'scripts/ecotype_pca_v2')); from fixed_projection_lib import stable_int_seed,VALID_GENOTYPES
  self.assertEqual(stable_int_seed('x',1),stable_int_seed('x',1)); self.assertEqual(VALID_GENOTYPES,frozenset('0129'))
