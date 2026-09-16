"""Include peptides whose sequence OR five-residue processing flanks contain17.

Unchanged peptide NetMHCpan scores come from the canonical full-protein scan;
changed peptides come from the20-choice screen. No new API calls are needed.
"""
import csv,json
from pathlib import Path
HERE=Path(__file__).resolve().parent;OUT=HERE/'mhcflurry_all_variants';NET=HERE.parent/'junction'
def read(p):return list(csv.DictReader(p.open(),delimiter='\t' if p.suffix=='.tsv' else ','))
def write(p,rows):
 with p.open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
base={}
variants={}
for m in ['ba','el']:
 base[m]={(r['allele'],int(r['start']),int(r['end'])):r for r in read(HERE/f'netmhcpan_{m}-4.1.tsv') if r['seq_num']=='2'}
 variants[m]={(r['allele'],int(r['seq_num']),int(r['start']),int(r['end'])):r for r in read(NET/f'variants_netmhcpan_{m}-4.1.tsv')}
substitutions=json.loads((NET/'variant_manifest.json').read_text())['construct_contexts']
seqs=json.loads((OUT/'input_sequences.json').read_text())
name_to_variant={name:'T17'+seq[16] for name,seq in seqs.items()}
name_to_num={name:next(i for i,x in enumerate(substitutions,1) if x['substitution']==v) for name,v in name_to_variant.items()}
allrows=read(OUT/'with_flanks.csv');rows=[]
for r in allrows:
 start=int(r['pos'])+1;end=start+len(r['peptide'])-1
 if not max(1,start-5)<=17<=end+5:continue
 a=r['sample_name'];n=name_to_num[r['sequence_name']];changed=start<=17<=end
 key=(a,n,start,end) if changed else (a,start,end)
 data=variants if changed else base
 ba,el=data['ba'][key],data['el'][key]
 assert ba['peptide']==el['peptide']==r['peptide']
 rows.append(dict(variant=name_to_variant[r['sequence_name']],allele=a,start=start,end=end,peptide=r['peptide'],
    peptide_contains17=changed,netmhcpan_el_percentile=float(el['percentile_rank']),netmhcpan_ba_nm=float(ba['ic50']),
    mhcflurry_ba_nm=float(r['affinity']),mhcflurry_presentation_percentile=float(r['presentation_percentile'])))
write(OUT/'peptide_or_flank_affected_windows.csv',rows)
metrics=[]
for v in dict.fromkeys(r['variant'] for r in rows):
 subset=[r for r in rows if r['variant']==v]
 result={'variant':v,'n':len(subset)}
 for t in [.5,1,2]:
  result[f'net_el_lt{t}']=sum(r['netmhcpan_el_percentile']<t for r in subset)
  result[f'mhc_presentation_lt{t}']=sum(r['mhcflurry_presentation_percentile']<t for r in subset)
 metrics.append(result)
keys=[f'{m}{t}' for m in ['net_el_lt','mhc_presentation_lt'] for t in [.5,1,2]]
def dominates(a,b):return all(a[k]<=b[k] for k in keys) and any(a[k]<b[k] for k in keys)
for r in metrics:
 r['dominated_by']=';'.join(s['variant'] for s in metrics if dominates(s,r))
 r['pareto']=not r['dominated_by']
write(OUT/'pareto_including_flank_affected_windows.csv',metrics)
print('Frontier',[r['variant'] for r in metrics if r['pareto']])
for r in metrics:
 if r['pareto'] or r['variant']=='T17T':print(r)
