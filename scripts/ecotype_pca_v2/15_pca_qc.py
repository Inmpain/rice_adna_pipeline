#!/usr/bin/env python3
import argparse
from fixed_projection_lib import read_evec,read_ind,write_tsv

p=argparse.ArgumentParser()
p.add_argument('--evec',required=True)
p.add_argument('--ind',required=True)
p.add_argument('--expected-n',type=int)
p.add_argument('--ancient-poplabel',default='Ancient',
               help='label used for appended ancient samples in --ind (matches '
                    '13_merge_ancients_fixed_panel.py --ancient-poplabel default). '
                    'smartpca lsqproject can silently drop a sample with near-zero '
                    'marker overlap from its own .evec output -- this is tolerated '
                    'ONLY for samples carrying this label, never for a modern '
                    'reference/axis-building sample.')
p.add_argument('--out',required=True)
a=p.parse_args()

e=read_evec(a.evec,10)
i=read_ind(a.ind)
evec_ids=[x['id'] for x in e]
ind_by_id={x['id']:x for x in i}

if len(evec_ids)!=len(set(evec_ids)):
 raise SystemExit('duplicate sample IDs in evec')

unknown_in_evec=sorted(set(evec_ids)-set(ind_by_id))
if unknown_in_evec:
 raise SystemExit(f'evec contains IDs not present in .ind: {unknown_in_evec[:10]}')

missing_from_evec=sorted(set(ind_by_id)-set(evec_ids))
missing_non_ancient=[s for s in missing_from_evec if ind_by_id[s]['label']!=a.ancient_poplabel]
if missing_non_ancient:
 raise SystemExit(f'{len(missing_non_ancient)} non-ancient sample(s) missing from evec '
                   f'(must never happen for a modern reference/axis sample): {missing_non_ancient[:10]}')

if a.expected_n is not None and len(i)!=a.expected_n:
 raise SystemExit(f'unexpected .ind sample count: {len(i)} != {a.expected_n}')

rows=[
 {'metric':'ind_sample_n','value':len(i)},
 {'metric':'evec_sample_n','value':len(evec_ids)},
 {'metric':'pc_n','value':10},
 {'metric':'dropped_by_smartpca_n','value':len(missing_from_evec)},
 {'metric':'dropped_by_smartpca_ids','value':','.join(missing_from_evec) if missing_from_evec else 'NONE'},
]
write_tsv(a.out,rows,['metric','value'])
print(f"PASS: evec={len(evec_ids)} ind={len(i)} dropped_by_smartpca={missing_from_evec or 'NONE'}")
