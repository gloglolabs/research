"""Double-mutant criterion audit including all processing-flank effects.

Positions9 and17 influence8-11mer peptide or five-residue flank contexts with
starts1-22. Net scores outside the altered sequence are reused from canonical
full protein only after exact peptide identity assertion.
"""
import csv,json,math
from pathlib import Path
HERE=Path(__file__).resolve().parent;OUT=HERE/'mhcflurry_leader_doubles_full';NET=HERE.parent/'leader_design/double_scan'
def read(p):return list(csv.DictReader(p.open(),delimiter='\t' if p.suffix=='.tsv' else ','))
def write(p,rows):
 with p.open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
man=json.loads((NET/'scan_manifest.json').read_text());constructs=man['constructs'];ids={c['id']:i for i,c in enumerate(constructs,1)}
seqs={c['id']:c['full_sequence'] for c in constructs};subs={c['id']:c['substitution'] for c in constructs}
net={m:{(int(r['seq_num']),r['allele'],int(r['start']),int(r['end'])):r for r in read(NET/f'netmhcpan_{m}-4.1.tsv')} for m in ['ba','el']}
base={m:{(r['allele'],int(r['start']),int(r['end'])):r for r in read(HERE/f'netmhcpan_{m}-4.1.tsv') if r['seq_num']=='2'} for m in ['ba','el']}
rows=[];checks={}
for context in ['with_flanks','without_flanks']:
 raw=read(OUT/f'{context}.csv');seen=set()
 for r in raw:
  start=int(r['pos'])+1;end=start+len(r['peptide'])-1
  if start>22:continue
  key=(ids[r['sequence_name']],r['sample_name'],start,end)
  scores={}
  for m in ['ba','el']:
   if key in net[m]:scores[m]=net[m][key]
   else:
    assert not (start<=9<=end or start<=17<=end)
    scores[m]=base[m][key[1:]]
   assert scores[m]['peptide']==r['peptide']
  assert key not in seen;seen.add(key)
  rows.append(dict(context=context,construct=r['sequence_name'],substitution=subs[r['sequence_name']],allele=r['sample_name'],start=start,end=end,peptide=r['peptide'],
    net_el_rank=float(scores['el']['percentile_rank']),mhc_presentation_rank=float(r['presentation_percentile']),mhc_affinity=float(r['affinity'])))
 assert len(seen)==len(constructs)*len(man['alleles'])*22*4
 checks[context+'_complete_extended_grid']=len(seen)
write(OUT/'processing_flanks_extended_joined.csv',rows)
metrics=[];totals=[]
for context in ['with_flanks','without_flanks']:
 vectors={}
 for c in constructs:
  vector={}
  for a in man['alleles']:
   selected=[r for r in rows if r['context']==context and r['construct']==c['id'] and r['allele']==a]
   for model in ['net_el_rank','mhc_presentation_rank']:
    for t in [.5,1,2]:vector[f'{a}|{model}|{t}']=sum(r[model]<t for r in selected)
  vectors[c['id']]=vector
  total={'context':context,'substitution':c['substitution'],'windows':616}
  for m in ['net_el_rank','mhc_presentation_rank']:
   for t in [.5,1,2]:total[f'{m}_lt{t}']=sum(vector[f'{a}|{m}|{t}'] for a in man['alleles'])
  totals.append(total)
 baseline=vectors[constructs[0]['id']]
 for c in constructs:
  vector=vectors[c['id']]
  metrics.append(dict(context=context,substitution=c['substitution'],**vector,better=sum(n<baseline[k] for k,n in vector.items()),same=sum(n==baseline[k] for k,n in vector.items()),worse=sum(n>baseline[k] for k,n in vector.items())))
write(OUT/'processing_flanks_extended_42_metrics.csv',metrics);write(OUT/'processing_flanks_extended_totals.csv',totals)
(OUT/'processing_flanks_extended_validation.json').write_text(json.dumps(checks,indent=2))
for r in metrics:print(r['context'],r['substitution'],r['better'],r['same'],r['worse'])
print('CC totals',[r for r in totals if r['substitution'] in ['baseline','L9C_T17C']])
