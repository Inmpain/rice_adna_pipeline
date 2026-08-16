#!/usr/bin/env python3
import argparse,collections
from fixed_projection_lib import read_ind,parse_sample_paths,read_calls,write_tsv
p=argparse.ArgumentParser(); p.add_argument('--ind',required=True); p.add_argument('--calls',action='append',required=True); p.add_argument('--folds',type=int,default=5); p.add_argument('--out',required=True); a=p.parse_args()
if a.folds!=5: raise SystemExit('exact-mask validation requires 5 folds')
ind=read_ind(a.ind); labels={x['id']:x['label'] for x in ind}; paths=parse_sample_paths(a.calls); counts=collections.Counter(labels.values()); rows=[]
for sample,path in paths:
 if sample not in labels: raise SystemExit(f'unknown sample {sample}')
 c=read_calls(path); rows.append({'sample':sample,'label':labels[sample],'fold':sum(map(ord,sample))%a.folds,'callable_n':sum(x!='9' for x in c),'status':'OK' if counts[labels[sample]]>=10 else 'INSUFFICIENT_REFERENCE_N'})
write_tsv(a.out,rows,['sample','label','fold','callable_n','status'])
