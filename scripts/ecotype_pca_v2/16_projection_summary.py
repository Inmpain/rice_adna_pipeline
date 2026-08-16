#!/usr/bin/env python3
import argparse
from fixed_projection_lib import read_evec,nearest_group,write_tsv

DIMENSIONS=(3,5,10)
K_VALUES=(20,50)

p=argparse.ArgumentParser(); p.add_argument('--panel',required=True,choices='ABC'); p.add_argument('--evec',required=True); p.add_argument('--out',required=True); a=p.parse_args()
rows=read_evec(a.evec,10); modern=[r for r in rows if not r['label'].lower().startswith('ancient')]; out=[]
for r in rows:
 d={'sample':r['id'],'label':r['label']}
 for dimensions in DIMENSIONS:
  for k in K_VALUES:
   try:
    winner,_max_count,_detail=nearest_group(modern,a.panel,r['pcs'],dimensions,k)
   except ValueError:
    winner='NA'
   d[f'nearest_pc{dimensions}_k{k}']=winner
 out.append(d)
fieldnames=['sample','label']+[f'nearest_pc{d}_k{k}' for d in DIMENSIONS for k in K_VALUES]
write_tsv(a.out,out,fieldnames)
