"""Recompute all plotted values from cached inputs; no network or model inference."""
from pathlib import Path
from collections import Counter, defaultdict
import argparse, csv, hashlib, json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
PRED = ROOT / 'novel/predictions'
DOUBLE = PRED / 'mhcflurry_leader_doubles_full'
NET = ROOT / 'novel/leader_design/double_scan'
FILES = {}
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--bundle-only', action='store_true', help='Ignore the original research tree; verify portable package inputs')
args = parser.parse_args()
BUNDLE = HERE.parent / 'reproducibility'
BUNDLE_MAP_PATH = BUNDLE / 'provenance/bundle_origins.json'
BUNDLE_MAP = {r['research_path']: (BUNDLE/r['bundle_path'], r['source_sha256']) for r in json.loads(BUNDLE_MAP_PATH.read_text())} if BUNDLE_MAP_PATH.exists() else {}
LOCAL_MAP = {r['research_path']: (HERE/r['figure_input_path'], r['source_sha256']) for r in json.loads((HERE/'source_inputs/origins.json').read_text())}

def resolve(path):
    key = path.relative_to(ROOT).as_posix()
    if not args.bundle_only and path.exists():
        resolved = path
    elif key in BUNDLE_MAP:
        resolved, expected_hash = BUNDLE_MAP[key]
        assert hashlib.sha256(resolved.read_bytes()).hexdigest() == expected_hash
    elif key in LOCAL_MAP:
        resolved, expected_hash = LOCAL_MAP[key]
        assert hashlib.sha256(resolved.read_bytes()).hexdigest() == expected_hash
    else:
        raise FileNotFoundError(f'No packaged source for {key}')
    FILES[key] = resolved
    return resolved

def read(path):
    resolved = resolve(path)
    return list(csv.DictReader(resolved.open(), delimiter='\t' if path.suffix == '.tsv' else ','))

def read_json(path):
    return json.loads(resolve(path).read_text())

def write(name, rows):
    with (HERE / name).open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)

ins = read_json(ROOT / 'data/uniprot_P01308.json')['sequence']['value']
chga = read_json(ROOT / 'data/uniprot_P10645.json')['sequence']['value']
man = read_json(NET / 'scan_manifest.json')
constructs = {c['substitution']: c for c in man['constructs']}
assert constructs['baseline']['full_sequence'] == chga[:18] + ins[24:]
selected = constructs['L9C_T17C']['full_sequence']
assert selected[:18] == 'MRSAAVLACLLCAGQVCA'
assert selected[18:] == ins[24:]
write('figure1_sequences.csv', [
    dict(construct='Native INS', sequence=ins, signal_end=24, proinsulin_source='INS25-110', type='native'),
    dict(construct='Canonical CHGA-INS', sequence=constructs['baseline']['full_sequence'], signal_end=18, proinsulin_source='INS25-110', type='reconstruction'),
    dict(construct='L9C + T17C', sequence=selected, signal_end=18, proinsulin_source='INS25-110', type='computed_candidate'),
])
lead = read_json(ROOT / 'novel/leader_design/scan_manifest.json')
t17 = read_json(ROOT / 'novel/junction/variant_manifest.json')
assert len(lead['constructs']) == 134 and len(t17['construct_contexts']) == 20
variant_names = [f'L9{x}_T17{y}' for x in ['C', 'G'] for y in ['C', 'G', 'D', 'P']]
assert all(v in constructs for v in variant_names)
write('figure1_scan_design.csv', [
    dict(stage='Junction site', positions='17', substitutions=19, controls=1, scored_sequences=20),
    dict(stage='N-terminal leader', positions='3-9', substitutions=133, controls=1, scored_sequences=134),
    dict(stage='Explicit combinations', positions='9 and 17', substitutions=8, controls=5, scored_sequences=13),
])

# Figure 2: counts independently recomputed from raw tool outputs.
alleles = ['HLA-A*02:01', 'HLA-A*03:01', 'HLA-A*24:02', 'HLA-B*40:01', 'HLA-B*49:01', 'HLA-C*03:04', 'HLA-C*07:01']
name_num = {i: c['substitution'] for i, c in enumerate(man['constructs'], 1)}
name_id = {c['id']: c['substitution'] for c in man['constructs']}
net = {}
for row in read(NET / 'netmhcpan_el-4.1.tsv'):
    if int(row['start']) > 18: continue
    key = (name_num[int(row['seq_num'])], row['allele'], int(row['start']), int(row['end']))
    assert key not in net
    net[key] = row
modes = {}
for mode in ['with_flanks', 'without_flanks']:
    lookup = {}
    for row in read(DOUBLE / f'{mode}.csv'):
        start = int(row['pos']) + 1; end = start + len(row['peptide']) - 1
        if start > 18: continue
        key = (name_id[row['sequence_name']], row['sample_name'], start, end)
        assert key not in lookup
        assert net[key]['peptide'] == row['peptide']
        lookup[key] = row
    assert lookup.keys() == net.keys()
    modes[mode] = lookup
assert len(net) == 13 * 7 * 72
summary_reference = {
    (r['substitution'], r['context'], r['allele'], float(r['threshold'])): r
    for r in read(DOUBLE / 'summary_counts.csv') if r['region'] == 'signal_or_junction'
}
counts = []
for variant in ['baseline'] + variant_names:
    for mode, lookup in modes.items():
        for allele in alleles:
            keys = [k for k in net if k[0] == variant and k[1] == allele]
            base_keys = [('baseline', k[1], k[2], k[3]) for k in keys]
            assert len(keys) == 72
            for cutoff in [.5, 1., 2.]:
                for model, field, source, ref in [
                    ('NetMHCpan EL', 'percentile_rank', net, 'net_el_lt'),
                    ('MHCflurry presentation', 'presentation_percentile', lookup, 'mhc_presentation_lt'),
                ]:
                    count = sum(float(source[k][field]) < cutoff for k in keys)
                    baseline_count = sum(float(source[k][field]) < cutoff for k in base_keys)
                    assert count == int(summary_reference[(variant, mode, allele, cutoff)][ref])
                    counts.append(dict(variant=variant, model=model, context=mode, allele=allele,
                        cutoff_percent=cutoff, windows=72, baseline_count=baseline_count,
                        variant_count=count, delta=count-baseline_count))
write('figure2_counts.csv', counts)
robustness = []
for variant in variant_names:
    for mode in modes:
        subset = [r for r in counts if r['variant'] == variant and r['context'] == mode]
        assert len(subset) == 42
        better=sum(r['delta']<0 for r in subset);same=sum(r['delta']==0 for r in subset);worse=sum(r['delta']>0 for r in subset)
        robustness.append(dict(variant=variant, context=mode, improved=better, unchanged=same, worsened=worse))
assert next(r for r in robustness if r['variant']=='L9C_T17C' and r['context']=='with_flanks')['improved']==16
assert next(r for r in robustness if r['variant']=='L9C_T17C' and r['context']=='without_flanks')['improved']==15
write('figure2_robustness.csv', robustness)
totals = []
for variant in ['baseline','L9C_T17C']:
    for model in ['NetMHCpan EL','MHCflurry presentation']:
        for mode in modes:
            for cutoff in [.5,1.,2.]:
                subset=[r for r in counts if r['variant']==variant and r['model']==model and r['context']==mode and r['cutoff_percent']==cutoff]
                total=sum(r['variant_count'] for r in subset)
                assert len(subset)==7
                totals.append(dict(variant=variant, model=model, context=mode, cutoff_percent=cutoff, windows=504, count=total))
write('figure2_totals.csv', totals)

# Figure 3: source observations and sequence coordinates are checked independently.
observations = [r for r in read(ROOT/'novel/peptidome/annotated_INS_CHGA_mhci.csv') if r['gene']=='INS']
donor_source = read(ROOT/'novel/peptidome/donor_footprint.csv')
peptides=['ALWGPDPAAA','AAAFVNQHL','GSHLVEALY','HLVEALYLV']
assert len(observations)==6 and set(r['peptide'] for r in observations)==set(peptides)
seen=set()
for r in observations:
    assert ins[int(r['start'])-1:int(r['end'])]==r['peptide']
    key=(r['peptide'],r['donor'].replace(' ',''));assert key not in seen;seen.add(key)
assert seen=={('ALWGPDPAAA','R361'),('ALWGPDPAAA','R369'),('AAAFVNQHL','HP20289'),('HLVEALYLV','HP20289'),('HLVEALYLV','HP18101'),('GSHLVEALY','R361')}
figure3=[]
for pep in peptides:
    original=next(r for r in observations if r['peptide']==pep)
    for d in donor_source:
        norm=d['donor'].replace(' ','');present=(pep,norm) in seen
        obs=next((r for r in observations if r['peptide']==pep and r['donor'].replace(' ','')==norm),None)
        figure3.append(dict(peptide=pep,start=int(original['start']),end=int(original['end']),
            footprint='retained' if int(original['start'])>24 else 'changed',donor=norm,detected=int(present),
            total_curated_classI_species=int(d['all_curated_mhci_species']),
            source_excel_row=obs['excel_row'] if obs else '',blank_intensity=obs['S4_exact_blank_intensity'] if obs else '',
            predicted_hla_as_published=obs['predicted_hla_as_published'] if obs else ''))
write('figure3_donor_projection.csv',figure3)

# Figure 4: quantitative cross-model checks, with the same HLA for each pair.
reference_ba = read(PRED/'netmhcpan_ba-4.1.tsv')
reference_mhc = read(PRED/'mhcflurry/with_flanks.csv')
double_ba = {}
for r in read(NET/'netmhcpan_ba-4.1.tsv'):
    key=(name_num[int(r['seq_num'])],r['allele'],int(r['start']),int(r['end']))
    double_ba[key]=r
double_rows=[]
for key,r in modes['with_flanks'].items():
    b=double_ba[key];e=net[key]
    assert b['peptide']==r['peptide']==e['peptide']
    double_rows.append(dict(substitution=key[0],allele=key[1],peptide=r['peptide'],context='with_flanks',
        region='signal' if key[3]<=18 else 'junction',net_ba_nm=b['ic50'],net_el_rank=e['percentile_rank'],
        mhc_ba_nm=r['affinity'],mhc_presentation_rank=r['presentation_percentile']))
junction=[]
for label,pep,cid in [('Native INS','AAAFVNQHL','WT_INS'),('Canonical fusion','VTAFVNQHL','CANONICAL_CHGA18_INS25_RECONSTRUCTION')]:
    seqnum=1 if cid=='WT_INS' else 2
    b=next(r for r in reference_ba if int(r['seq_num'])==seqnum and r['peptide']==pep and r['allele']=='HLA-C*03:04')
    m=next(r for r in reference_mhc if r['sequence_name']==cid and r['peptide']==pep and r['sample_name']=='HLA-C*03:04')
    junction.append(dict(label=label,peptide=pep,allele=b['allele'],net_affinity_nm=float(b['ic50']),mhc_affinity_nm=float(m['affinity'])))
for label,pep in [('L9C + T17C','VCAFVNQHL'),('Nested candidate','CAFVNQHL')]:
    r=next(r for r in double_rows if r['substitution']=='L9C_T17C' and r['peptide']==pep and r['allele']=='HLA-C*03:04' and r['context']=='with_flanks')
    junction.append(dict(label=label,peptide=pep,allele=r['allele'],net_affinity_nm=float(r['net_ba_nm']),mhc_affinity_nm=float(r['mhc_ba_nm'])))
write('figure4_junction_affinity.csv',junction)
residual=[]
for r in double_rows:
    if r['substitution']=='L9C_T17C' and r['context']=='with_flanks' and float(r['mhc_presentation_rank'])<.5:
        residual.append(dict(peptide=r['peptide'],allele=r['allele'],region=r['region'],net_el_rank_percent=float(r['net_el_rank']),mhc_presentation_rank_percent=float(r['mhc_presentation_rank']),net_affinity_nm=float(r['net_ba_nm']),mhc_affinity_nm=float(r['mhc_ba_nm'])))
assert len(residual)==3 and {r['allele'] for r in residual}=={'HLA-C*03:04'}
write('figure4_residual_candidates.csv',residual)
trade=[]
for variant,pep in [('baseline','VLALLLCAGQV'),('L9C_T17C','VLACLLCAGQV'),('L9G_T17C','VLAGLLCAGQV')]:
    r=next(r for r in double_rows if r['substitution']==variant and r['peptide']==pep and r['allele']=='HLA-A*02:01' and r['context']=='with_flanks')
    trade.append(dict(variant=variant,peptide=pep,allele=r['allele'],net_affinity_nm=float(r['net_ba_nm']),mhc_affinity_nm=float(r['mhc_ba_nm']),net_el_rank_percent=float(r['net_el_rank']),mhc_presentation_rank_percent=float(r['mhc_presentation_rank'])))
write('figure4_A02_affinity_tradeoff.csv',trade)
manifest={
    'scope':'Original plots of cached predictions and public donor observations; no simulated data or new model runs',
    'source_files':[{'path':key,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()} for key,path in sorted(FILES.items())],
    'validation':{'net_and_mhc_complete_matching_signal_junction_grid':len(net),'all_recomputed_counts_match_saved_reference':True,'donor_observations':6,'distinct_INS_peptides':4,'INS_positive_donors':4,'donors_with_retained_INS':3,'all_INS_coordinates_match_UniProt':True},
    'figure_source_csvs':[p.name for p in sorted(HERE.glob('figure*.csv'))],
}
(HERE/'source_data_manifest.json').write_text(json.dumps(manifest,indent=2))
print(json.dumps(manifest['validation'],indent=2))
