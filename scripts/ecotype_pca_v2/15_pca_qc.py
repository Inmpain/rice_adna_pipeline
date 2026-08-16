#!/usr/bin/env python3
import argparse
from fixed_projection_lib import read_evec,read_ind,write_tsv
p=argparse.ArgumentParser(); p.add_argument('--evec',required=True); p.add_argument('--ind',required=True); p.add_argument('--expected-n',type=int); p.add_argument('--out',required=True); a=p.parse_args()
e=read_evec(a.evec,10); i=read_ind(a.ind); ids=[x['id'] for x in e]; ind_ids=[x['id'] for x in i]
if len(ids)!=len(set(ids)): raise SystemExit('duplicate IDs in evec')
unknown=sorted(set(ids)-set(ind_ids))
if unknown: raise SystemExit('evec contains IDs absent from ind: '+','.join(unknown[:10]))
missing=sorted(set(ind_ids)-set(ids))
if a.expected_n is not None and len(ids)>a.expected_n: raise SystemExit('evec has more samples than expected')
write_tsv(a.out,[{'metric':'evec_sample_n','value':len(ids)},{'metric':'ind_sample_n','value':len(i)},{'metric':'evec_samples_missing_from_pca','value':len(missing)},{'metric':'missing_sample_ids','value':','.join(missing)},{'metric':'pc_n','value':10}],['metric','value'])
