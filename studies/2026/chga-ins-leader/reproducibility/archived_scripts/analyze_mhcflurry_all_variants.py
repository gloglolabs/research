"""Join root's NetMHCpan T17 screen to independent full-protein MHCflurry checks."""
import csv,json
from pathlib import Path
HERE=Path(__file__).resolve().parent;OUT=HERE/'mhcflurry_all_variants';NET=HERE.parent/'junction'
def read(p):return list(csv.DictReader(p.open(),delimiter='\t' if p.suffix=='.tsv' else ','))
def write(p,rows):
 with p.open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
variants=json.loads((NET/'variant_manifest.json').read_text())['construct_contexts']
methods={m:{(r['allele'],int(r['seq_num']),int(r['start']),int(r['end'])):r for r in read(NET/f'variants_netmhcpan_{m}-4.1.tsv')} for m in ['ba','el']}
constructs=json.loads((OUT/'input_sequences.json').read_text())
lookup={}
for name,seq in constructs.items():
 x=seq[16];n=next(i for i,v in enumerate(variants,1) if v['substitution']==f'T17{x}')
 assert seq.startswith(variants[n-1]['sequence'])
 lookup[name]=n
rows=[];validation={}
for mode in ['with_flanks','without_flanks']:
 raw=read(OUT/f'{mode}.csv');seen=set()
 for r in raw:
  start=int(r['pos'])+1;end=start+len(r['peptide'])-1
  key=(r['sequence_name'],r['sample_name'],start,end)
  assert key not in seen;seen.add(key)
  assert constructs[r['sequence_name']][start-1:end]==r['peptide']
  if not start<=17<=end:continue
  nk=(r['sample_name'],lookup[r['sequence_name']],start,end);ba=methods['ba'][nk];el=methods['el'][nk]
  assert ba['peptide']==el['peptide']==r['peptide']
  rows.append(dict(context=mode,variant=variants[nk[1]-1]['substitution'],allele=r['sample_name'],start=start,end=end,
    peptide=r['peptide'],netmhcpan_ba_nm=float(ba['ic50']),netmhcpan_ba_percentile=float(ba['percentile_rank']),
    netmhcpan_el_percentile=float(el['percentile_rank']),mhcflurry_ba_nm=float(r['affinity']),
    mhcflurry_ba_percentile=float(r['affinity_percentile']),mhcflurry_processing_score=float(r['processing_score']),
    mhcflurry_presentation_percentile=float(r['presentation_percentile'])))
 expected={(name,a,start,start+n-1) for name,seq in constructs.items() for a in json.loads((OUT/'manifest.json').read_text())['alleles'] for n in [8,9,10,11] for start in range(1,len(seq)-n+2)}
 assert seen==expected
 validation[mode+'_complete_grid_rows']=len(seen)
write(OUT/'affected_windows_joined.csv',rows)
summary=[]
for mode in ['with_flanks','without_flanks']:
 for v in [f'T17{x}' for x in 'ACDEFGHIKLMNPQRSTVWY']:
  for a in ['all']+list(json.loads((OUT/'manifest.json').read_text())['alleles']):
   selected=[r for r in rows if r['context']==mode and r['variant']==v and (a=='all' or r['allele']==a)]
   for t in [.5,1,2,5]:
    summary.append(dict(context=mode,variant=v,allele=a,threshold=t,n_windows=len(selected),
       netmhcpan_el_lt=sum(r['netmhcpan_el_percentile']<t for r in selected),
       mhcflurry_presentation_lt=sum(r['mhcflurry_presentation_percentile']<t for r in selected),
       both_presentation_lt=sum(r['netmhcpan_el_percentile']<t and r['mhcflurry_presentation_percentile']<t for r in selected)))
write(OUT/'threshold_summary.csv',summary)
write(OUT/'highlight_peptides.csv',[r for r in rows if r['start']==16 and r['end']==24])
validation['affected_windows_joined_per_context']=len(rows)//2
(OUT/'validation.json').write_text(json.dumps(validation,indent=2));print(json.dumps(validation))
