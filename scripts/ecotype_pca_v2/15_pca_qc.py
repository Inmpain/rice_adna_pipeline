#!/usr/bin/env python3
import argparse
from fixed_projection_lib import read_evec,read_ind,write_tsv
p=argparse.ArgumentParser(); p.add_argument('--evec',required=True); p.add_argument('--ind',required=True); p.add_argument('--expected-n',type=int); p.add_argument('--out',required=True); a=p.parse_args()
e=read_evec(a.evec,10); i=read_ind(a.ind); ids=[x['id'] for x in e]
if len(ids)!=len(set(ids)) or len(ids)!=len(i): raise SystemExit('sample count/ID mismatch')
if a.expected_n is not None and len(ids)!=a.expected_n: raise SystemExit('unexpected sample count')
write_tsv(a.out,[{'metric':'sample_n','value':len(ids)},{'metric':'pc_n','value':10}],['metric','value'])
