#!/usr/bin/env python3
import argparse
from fixed_projection_lib import read_evec,class_centroids,nearest_group,write_tsv
p=argparse.ArgumentParser(); p.add_argument('--panel',required=True,choices='ABC'); p.add_argument('--evec',required=True); p.add_argument('--out',required=True); a=p.parse_args()
rows=read_evec(a.evec,10); modern=[r for r in rows if not r['label'].lower().startswith('ancient')]; out=[]
for r in rows:
 d={'sample':r['id'],'label':r['label']}
 for n in (3,5,10): d[f'nearest_pc{n}']=nearest_group(r['pcs'],class_centroids(modern,a.panel,n),n) or 'NA'
 out.append(d)
write_tsv(a.out,out,['sample','label','nearest_pc3','nearest_pc5','nearest_pc10'])
