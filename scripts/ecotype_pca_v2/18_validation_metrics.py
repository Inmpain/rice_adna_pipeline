#!/usr/bin/env python3
import argparse,csv
from fixed_projection_lib import write_tsv
p=argparse.ArgumentParser(); p.add_argument('--assignments',required=True); p.add_argument('--min-class-n',type=int,default=10); p.add_argument('--out',required=True); a=p.parse_args(); rows=list(csv.DictReader(open(a.assignments),delimiter='\t'))
truth=lambda r:r.get('truth',r.get('label','')); pred=lambda r:r.get('predicted',r.get('assignment','')); out=[]
for c in sorted(set(truth(r) for r in rows)):
 tp=sum(truth(r)==c and pred(r)==c for r in rows); n=sum(truth(r)==c for r in rows); pp=sum(pred(r)==c for r in rows)
 out.append({'class':c,'n':n,'precision':tp/pp if pp else 0.0,'recall':tp/n if n else 0.0,'status':'OK' if n>=a.min_class_n else 'INSUFFICIENT_REFERENCE_N'})
write_tsv(a.out,out,['class','n','precision','recall','status'])
