#!/usr/bin/env python3
"""Reproduce manuscript tables from cached predictions and attributed source excerpts.
Python>=3.10, standard library only. No network, model weights, or project checkout.
"""
import argparse,csv,hashlib,json,math
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parent
ALLELES=['HLA-A*02:01','HLA-A*24:02','HLA-A*03:01','HLA-B*40:01','HLA-B*49:01','HLA-C*03:04','HLA-C*07:01']
LENGTHS=[8,9,10,11]
CUTS=[0.5,1,2]

def read_csv(p,delimiter=','):
    with p.open(newline='') as f:return list(csv.DictReader(f,delimiter=delimiter))

def json_read(p):return json.loads(p.read_text())

def csv_write(p,rows,fields=None):
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields or list(rows[0]));w.writeheader();w.writerows(rows)

def fasta_read(path):
    result={};name=None
    for line in path.read_text().splitlines():
        if line.startswith('>'):name=line[1:].split()[0];assert name not in result;result[name]=''
        elif line.strip():result[name]+=line.strip()
    return result

def sha256(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()

def check_integrity():
    manifest=ROOT/'checksums.sha256'
    assert manifest.is_file(),'Missing checksums.sha256'
    count=0
    for line in manifest.read_text().splitlines():
        digest,name=line.split('  ',1)
        assert sha256(ROOT/name)==digest,f'Checksum mismatch: {name}'
        count+=1
    return count

def stage_analysis(stage,out):
    m=json_read(ROOT/f'metadata/{stage}_scan.json')
    cs=m.get('constructs',m.get('construct_contexts'))
    seqs={i+1:c['sequence'] for i,c in enumerate(cs)}
    ids={i+1:c['id'] for i,c in enumerate(cs)}
    d={mode:{i:{} for i in seqs} for mode in ['el','ba']}
    for mode in d:
        for r in read_csv(ROOT/f'cached/{stage}/netmhcpan_{mode}.tsv','\t'):
            i=int(r['seq_num']);s=int(r['start']);e=int(r['end']);a=r['allele'];k=(a,s,e)
            assert a in ALLELES and e-s+1 in LENGTHS
            assert seqs[i][s-1:e]==r['peptide'] and len(r['peptide'])==int(r['length'])
            assert k not in d[mode][i]
            d[mode][i][k]=r
        for i,seq in seqs.items():
            expected={(a,s,s+n-1) for a in ALLELES for n in LENGTHS for s in range(1,len(seq)-n+2)}
            assert set(d[mode][i])==expected
    rows=[]
    for i,c in enumerate(cs,1):
        er=d['el'][i];br=d['ba'][i]
        variant=c.get('substitution',c['id'].replace('CHGA_',''))
        for scope,pred in [('leader',lambda s,e:e<=18),('junction',lambda s,e:s<=18<e),('signal_or_junction',lambda s,e:s<=18),('position17_containing',lambda s,e:s<=17<=e)]:
            keys=[k for k in er if pred(k[1],k[2])]
            r={'stage':stage,'construct_id':ids[i],'substitution':variant,'scope':scope,'window_HLA_pairs':len(keys)}
            for cut in CUTS:r[f'EL_lt_{cut}']=sum(float(er[k]['percentile_rank'])<cut for k in keys)
            r['BA_lt_500nM']=sum(float(br[k]['ic50'])<500 for k in keys)
            rows.append(r)
    csv_write(out/f'{stage}_netmhcpan_summary.csv',rows)
    return cs,d,rows



def reference_comparison(stage_results, out):
    """Rebuild the three selected junction peptide comparisons from raw scores."""
    constructs, net, _ = stage_results["reference"]
    sequences = json_read(ROOT / "metadata/reference_mhcflurry_sequences.json")
    data = {}
    for mode in ["with_flanks", "without_flanks"]:
        data[mode] = {}
        for row in read_csv(ROOT / f"cached/reference/mhcflurry_{mode}.csv"):
            name = row["sequence_name"]
            start = int(row["pos"]) + 1
            end = start + len(row["peptide"]) - 1
            key = (name, row["sample_name"], start, end)
            assert key not in data[mode]
            assert sequences[name][start-1:end] == row["peptide"]
            data[mode][key] = row
        expected = {(name, a, s, s+n-1) for name, seq in sequences.items()
                    for a in ALLELES for n in LENGTHS for s in range(1, len(seq)-n+2)}
        assert set(data[mode]) == expected
    table = []
    targets = {"WT_INS": ["AAAFVNQHL"],
               "CANONICAL_CHGA18_INS25_RECONSTRUCTION": ["VTAFVNQHL", "TAFVNQHL"]}
    for index, construct in enumerate(constructs, 1):
        name = construct["id"]
        for peptide in targets.get(name, []):
            start = sequences[name].index(peptide) + 1
            end = start + len(peptide) - 1
            for allele in ALLELES:
                key = (allele, start, end)
                el = net["el"][index][key]
                ba = net["ba"][index][key]
                mfkey = (name, allele, start, end)
                with_flanks = data["with_flanks"][mfkey]
                without_flanks = data["without_flanks"][mfkey]
                table.append({"construct": name, "peptide": peptide, "HLA": allele,
                              "start_1based": start, "end_1based": end,
                              "net_EL_percentile": el["percentile_rank"],
                              "net_BA_nM": ba["ic50"],
                              "mhc_affinity_nM": with_flanks["affinity"],
                              "mhc_presentation_percentile_with_flanks": with_flanks["presentation_percentile"],
                              "mhc_presentation_percentile_without_flanks": without_flanks["presentation_percentile"]})
    assert len(table) == 21
    values = {r["peptide"]: r for r in table if r["HLA"] == "HLA-C*03:04"}
    assert float(values["AAAFVNQHL"]["net_BA_nM"]) == 46.38
    assert float(values["VTAFVNQHL"]["net_BA_nM"]) == 263.22
    assert float(values["TAFVNQHL"]["net_EL_percentile"]) == 0.27
    csv_write(out / "reference_junction_comparison.csv", table)
    return {"reference_MHCflurry_sequences": len(sequences), "selected_comparison_rows": len(table)}

def single_mhc_analysis(stage_results, out, refs):
    """Validate both single scans; reconstruct the position-17 sensitivity grid."""
    all_summaries = []
    stage_counts = {}
    position17_raw = {}
    position17_sequences = {}
    for stage in ["leader_3_9", "position17"]:
        constructs, _, _ = stage_results[stage]
        if stage == "leader_3_9":
            sequences = {c["id"]: c["full_sequence"] for c in constructs}
        else:
            sequences = {
                "CANONICAL_CHGA18_" + c["substitution"] + "_INS25_SENSITIVITY":
                c["sequence"][:18] + refs["INS"]["sequence"][24:]
                for c in constructs
            }
            position17_sequences = sequences
        for mode in ["with_flanks", "without_flanks"]:
            data = {name: {} for name in sequences}
            rows = read_csv(ROOT / f"cached/{stage}/mhcflurry_{mode}.csv")
            stage_counts[f"{stage}_{mode}"] = len(rows)
            for row in rows:
                name = row["sequence_name"]
                sequence = sequences[name]
                start = int(row["pos"]) + 1
                end = start + len(row["peptide"]) - 1
                key = (row["sample_name"], start, end)
                assert key not in data[name]
                assert sequence[start-1:end] == row["peptide"]
                if mode == "with_flanks":
                    assert row["n_flank"] == sequence[max(0, start-6):start-1]
                    assert row["c_flank"] == sequence[end:end+5]
                data[name][key] = row
            for name, sequence in sequences.items():
                expected = {(a, s, s+n-1) for a in ALLELES for n in LENGTHS
                            for s in range(1, len(sequence)-n+2)
                            if stage != "leader_3_9" or s <= 18}
                assert set(data[name]) == expected
                for scope, predicate in [
                    ("leader", lambda s, e: e <= 18),
                    ("signal_or_junction", lambda s, e: s <= 18),
                    ("position17_containing", lambda s, e: s <= 17 <= e),
                ]:
                    selected = [row for (_, s, e), row in data[name].items()
                                if predicate(s, e)]
                    record = {"stage": stage, "construct": name, "context": mode,
                              "scope": scope, "window_HLA_pairs": len(selected)}
                    for cut in CUTS:
                        record[f"presentation_lt_{cut}"] = sum(
                            float(row["presentation_percentile"]) < cut for row in selected)
                    all_summaries.append(record)
            if stage == "position17":
                position17_raw[mode] = data
    assert stage_counts["leader_3_9_with_flanks"] == 67536
    assert stage_counts["position17_with_flanks"] == 53480
    csv_write(out / "single_mhcflurry_summary.csv", all_summaries)

    # Peptide-only NetMHC scores beyond each 29-aa request can be reused only
    # after exact identity to the independently scanned full canonical precursor.
    constructs, net, _ = stage_results["position17"]
    reference_constructs, reference_net, _ = stage_results["reference"]
    reference_index = next(i for i, c in enumerate(reference_constructs, 1)
        if c["id"] == "CANONICAL_CHGA18_INS25_RECONSTRUCTION")
    name_index = {"CANONICAL_CHGA18_" + c["substitution"] + "_INS25_SENSITIVITY": i
                  for i, c in enumerate(constructs, 1)}
    baseline = "CANONICAL_CHGA18_T17T_INS25_SENSITIVITY"
    for name, index in name_index.items():
        for key, row in reference_net["el"][reference_index].items():
            if key not in net["el"][index]:
                assert not key[1] <= 17 <= key[2]
                assert position17_sequences[name][key[1]-1:key[2]] == row["peptide"]
                net["el"][index][key] = row
    objectives = []
    for mode in ["with_flanks", "without_flanks"]:
        for scope, predicate, denominator in [
            ("position17_containing", lambda s, e: s <= 17 <= e, 266),
            ("position17_core_or_flank", lambda s, e: s-5 <= 17 <= e+5, 546),
        ]:
            for name, index in name_index.items():
                for model in ["NetMHC", "MHCflurry"]:
                    data = net["el"][index] if model == "NetMHC" else position17_raw[mode][name]
                    base_data = net["el"][name_index[baseline]] if model == "NetMHC" else position17_raw[mode][baseline]
                    column = "percentile_rank" if model == "NetMHC" else "presentation_percentile"
                    selected = [key for key in data if predicate(key[1], key[2])]
                    assert len(selected) == denominator
                    for allele in ALLELES:
                        for cut in CUTS:
                            count = sum(float(data[key][column]) < cut for key in selected if key[0] == allele)
                            base_count = sum(float(base_data[key][column]) < cut for key in selected if key[0] == allele)
                            objectives.append({"variant": constructs[index-1]["substitution"],
                                "context": mode, "scope": scope, "window_HLA_pairs": denominator,
                                "model": model, "HLA": allele, "cutoff_percentile": cut,
                                "baseline_count": base_count, "variant_count": count,
                                "delta": count-base_count})
    summaries = []
    for mode in ["with_flanks", "without_flanks"]:
        for scope in ["position17_containing", "position17_core_or_flank"]:
            for c in constructs:
                selected = [r for r in objectives if r["variant"] == c["substitution"]
                            and r["context"] == mode and r["scope"] == scope]
                assert len(selected) == 42
                summaries.append({"variant": c["substitution"], "context": mode, "scope": scope,
                                  "better": sum(r["delta"] < 0 for r in selected),
                                  "unchanged": sum(r["delta"] == 0 for r in selected),
                                  "worse": sum(r["delta"] > 0 for r in selected)})
    feasible = {}
    for scope in ["position17_containing", "position17_core_or_flank"]:
        feasible[scope] = [r["variant"] for r in summaries if r["context"] == "with_flanks"
                           and r["scope"] == scope and r["better"] > 0 and r["worse"] == 0]
        assert feasible[scope] == ["T17C"]
    originals = {r["variant"]: r for r in read_csv(ROOT / "tables/position17_flank_affected_42_objective_metrics.csv")}
    for r in objectives:
        if r["scope"] == "position17_core_or_flank" and r["context"] == "with_flanks":
            original_model = "netmhcpan_el_percentile" if r["model"] == "NetMHC" else "mhcflurry_presentation_percentile"
            column = f"{r['HLA']}|{original_model}|{r['cutoff_percentile']}"
            assert r["variant_count"] == int(originals[r["variant"]][column])
    csv_write(out / "position17_objectives42.csv", objectives)
    csv_write(out / "position17_criteria42_summary.csv", summaries)
    first_signalp = fasta_read(ROOT / "cached/signalp_position17/processed_entries.fasta")
    assert len(first_signalp) == 21 and all(s == refs["INS"]["sequence"][24:] for s in first_signalp.values())
    return {"validated_raw_rows": stage_counts, "position17_core_pairs": 266,
            "position17_core_or_flank_pairs": 546, "only_improving_no_worse_with_flanks": feasible,
            "first_SignalP_batch_identical_mature_products": 21}


def inflammation_analysis(out, refs):
    """Recalculate construct projection and source-summary checks from excerpts."""
    rows = read_csv(ROOT / "tables/inflammation_ins_canonical_projection.csv")
    assert len(rows) == 29 and len({r["peptide"] for r in rows}) == 29
    groups = {}
    for row in rows:
        start, end = int(row["start"]), int(row["end"])
        assert refs["INS"]["sequence"][start-1:end] == row["peptide"]
        region = "INS signal 1-24" if end <= 24 else "INS boundary overlap" if start <= 24 else "INS retained 25-110"
        assert row["construct_region"] == region
        groups.setdefault(region, []).append(row)
        if row["mean_basal_intensity"] and row["mean_ifna_intensity"]:
            ratio = float(row["mean_ifna_intensity"]) / float(row["mean_basal_intensity"])
            assert math.isclose(ratio, float(row["ifna_fold_change"]), rel_tol=1e-12)
    totals = []
    for region, selected in sorted(groups.items()):
        totals.append({"region": region, "distinct_peptides": len(selected),
                       "basal_identified_at_least_2of4": sum(int(r["basal_identifications"]) >= 2 for r in selected),
                       "IFNa_identified_at_least_2of4": sum(int(r["ifna_identifications"]) >= 2 for r in selected)})
    lookup = {r["region"]: r for r in totals}
    assert lookup["INS retained 25-110"]["distinct_peptides"] == 25
    assert lookup["INS retained 25-110"]["basal_identified_at_least_2of4"] == 8
    assert lookup["INS retained 25-110"]["IFNa_identified_at_least_2of4"] == 12
    assert lookup["INS signal 1-24"]["distinct_peptides"] == 3
    assert lookup["INS boundary overlap"]["distinct_peptides"] == 1
    csv_write(out / "inflammation_projection_summary.csv", totals)
    run_rows = read_csv(ROOT / "tables/inflammation_figure8c_reported_run_intensities.csv")
    source_checks = read_csv(ROOT / "tables/inflammation_source_summary_consistency.csv")
    for check in source_checks:
        selected = [r for r in run_rows if r["peptide"] == check["peptide"]]
        assert len(selected) == 4
        for condition in ["basal", "ifna"]:
            mean = sum(float(r[f"{condition}_intensity"]) for r in selected) / len(selected)
            assert math.isclose(mean, float(check[f"figure8c_mean_{condition}"]), rel_tol=1e-12)
    return totals

def donor_analysis(out,refs):
    obs=read_csv(ROOT/'tables/donor_observations.csv')
    depths=read_csv(ROOT/'tables/donor_depth_and_projection.csv')
    ins=[r for r in obs if r['gene']=='INS'];graft=refs['CHGA']['sequence'][:18]+refs['INS']['sequence'][24:]
    assert len({(r['donor'],r['peptide']) for r in obs})==len(obs)==9
    for r in obs:
        s=int(r['start']);e=int(r['end']);seq=refs[r['gene']]['sequence']
        assert seq[s-1:e]==r['peptide'] and len(r['peptide'])==int(r['length'])
        assert (r['peptide'] in graft)==(r['sequence_in_canonical_CHGA1_18_INS25_110']=='True')
        if r['gene']=='INS':
            footprint='inside_native_INS_signal' if e<=24 else ('crosses_native_INS_signal_boundary' if s<=24 else 'downstream_INS_retained_sequence')
            assert footprint==r['footprint']
    assert len(ins)==6 and len({r['peptide'] for r in ins})==4
    donors={r['donor'] for r in depths};assert len(donors)==6
    assert sum(int(r['all_curated_mhci_species']) for r in depths)==3766
    def summarize(label,rows):
        kept=[r for r in rows if int(r['start'])>=25]
        return {'analysis':label,'observations':len(rows),'distinct_peptides':len({r['peptide'] for r in rows}),'INS_positive_donors':len({r['donor'] for r in rows}),'retained_positive_donors':len({r['donor'] for r in kept}),'retained_distinct_peptides':len({r['peptide'] for r in kept})}
    cases=[('authors_curated_S5',ins),('drop_exact_blank_positive_keep_unknown',[r for r in ins if not r['S4_exact_blank_intensity'] or float(r['S4_exact_blank_intensity'])==0]),('require_explicit_zero_exact_blank',[r for r in ins if r['S4_exact_blank_intensity']!='' and float(r['S4_exact_blank_intensity'])==0]),('drop_same_donor_6aa_blank_overlap_keep_unknown',[r for r in ins if not r['same_donor_blank_overlap_at_least_6aa']]),('exclude_lowest_depth_donor_HP18101_16_total_peptides',[r for r in ins if r['donor']!='HP18101'])]
    sensitivity=[summarize(label,rows) for label,rows in cases]
    assert sensitivity[0]['INS_positive_donors']==4 and sensitivity[0]['retained_positive_donors']==3
    reference=read_csv(ROOT/'tables/donor_background_sensitivity_original.csv')
    for a,b in zip(sensitivity,reference):
        assert a['analysis']==b['analysis']
        assert a['observations']==int(b['INS_donor_peptide_observations'])
        assert a['retained_positive_donors']==int(b['donors_with_retained_INS_species'])
    csv_write(out/'donor_sensitivity.csv',sensitivity)
    csv_write(out/'donor_leave_one_out.csv',[summarize('exclude_'+d,[r for r in ins if r['donor']!=d]) for d in sorted(donors)])
    peptides=sorted({r['peptide'] for r in ins},key=lambda p:refs['INS']['sequence'].index(p))
    matrix=[]
    for donor in sorted(donors):
        for p in peptides:
            s=refs['INS']['sequence'].index(p)+1;e=s+len(p)-1
            matched=[r for r in ins if r['donor']==donor and r['peptide']==p]
            matrix.append({'donor_as_published':donor,'donor_label':donor.replace(' ',''),'peptide':p,'INS_start_1based':s,'INS_end_1based':e,'detected_in_curated_table':bool(matched),'retained_in_canonical_graft':s>=25,'ion_intensity':matched[0]['ion_intensity'] if matched else '', 'predicted_HLA_as_published':matched[0]['predicted_hla_as_published'] if matched else ''})
    assert len(matrix) == 24 and sum(r['detected_in_curated_table'] for r in matrix) == 6
    csv_write(out/'donor_peptide_matrix.csv',matrix)
    return sensitivity[0]

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',default='results',help='Output directory, relative to bundle unless absolute');args=ap.parse_args()
    out=Path(args.output);out=out if out.is_absolute() else ROOT/out;out.mkdir(parents=True,exist_ok=True)
    file_count=check_integrity()
    refs=json_read(ROOT/'metadata/reference_proteins.json')
    stage_results={}
    for stage in ['leader_3_9','position17','double','reference']:
        cs,d,rows=stage_analysis(stage,out);stage_results[stage]=(cs,d,rows)
    assert len(stage_results['leader_3_9'][0])==134
    assert len(stage_results['position17'][0])==20
    zero_leader=[r['substitution'] for r in stage_results['leader_3_9'][2] if r['scope']=='leader' and r['EL_lt_0.5']==0]
    assert sorted(zero_leader)==['L9C','L9G']
    constructs,net,net_summary=stage_results['double'];assert len(constructs)==13
    cs={c['id']:c for c in constructs};ids=list(cs);base='CHGA_baseline';candidate='CHGA_L9C_T17C'
    sequences=fasta_read(ROOT/'constructs.fasta');assert sequences=={c['id']:c['full_sequence'] for c in constructs}
    mature=refs['INS']['sequence'][24:]
    for c in constructs:assert c['full_sequence'][18:]==mature and len(c['full_sequence'])==104
    assert [i+1 for i,(a,b) in enumerate(zip(sequences[base],sequences[candidate])) if a!=b]==[9,17]
    signalp=fasta_read(ROOT/'cached/signalp/input.fasta');processed=fasta_read(ROOT/'cached/signalp/processed_entries.fasta')
    assert all(signalp[k.upper()]==v for k,v in sequences.items())
    assert len(processed)==16 and all(s==mature for s in processed.values())
    sp=json_read(ROOT/'cached/signalp/validation.json')['predictions'];sp={r['id']:r for r in sp}
    assert sp['CHGA_L9C_T17C']['cleavage_after']==sp['CHGA_BASELINE']['cleavage_after']==18
    # Reconstruct complete NetMHC grids using the separately cached full canonical
    # precursor only for exactly identical peptide sequences beyond29aa.
    reference_cs,reference_net,_=stage_results['reference']
    ri=next(i for i,c in enumerate(reference_cs,1) if c['id']=='CANONICAL_CHGA18_INS25_RECONSTRUCTION')
    direct_net_keys={i:set(v) for i,v in net['el'].items()}
    reused_net_pairs=0
    for i,c in enumerate(constructs,1):
        for model in ['el','ba']:
            for k,rr in reference_net[model][ri].items():
                if k in net[model][i]:continue
                assert not any(k[1]<=p<=k[2] for p in c['positions'])
                assert c['full_sequence'][k[1]-1:k[2]]==rr['peptide']
                net[model][i][k]=rr
                if model=='el':reused_net_pairs+=1
            assert len(net[model][i])==2674
    assert reused_net_pairs==13*(2674-574)
    raw={'NetMHC':{ids[i-1]:v for i,v in net['el'].items()}}
    for mode in ['with_flanks','without_flanks']:
        raw[mode]={c:{} for c in ids}
        for r in read_csv(ROOT/f'cached/double/mhcflurry_{mode}.csv'):
            c=r['sequence_name'];a=r['sample_name'];s=int(r['pos'])+1;e=s+len(r['peptide'])-1;k=(a,s,e)
            assert a in ALLELES and e-s+1 in LENGTHS and sequences[c][s-1:e]==r['peptide']
            assert k not in raw[mode][c];raw[mode][c][k]=r
            if mode=='with_flanks':
                assert r['n_flank']==sequences[c][max(0,s-6):s-1] and r['c_flank']==sequences[c][e:e+5]
        for c in ids:
            expected={(a,s,s+n-1) for a in ALLELES for n in LENGTHS for s in range(1,105-n+1)}
            assert len(expected)==2674 and set(raw[mode][c])==expected
    # Outside core and relevant5-aa flanks, mutation cannot affect these inputs.
    for mode in ['with_flanks','without_flanks']:
        for k,w in raw[mode][candidate].items():
            radius=5 if mode=='with_flanks' else 0
            affected=any(k[1]-radius<=p<=k[2]+radius for p in [9,17])
            if not affected:
                b=raw[mode][base][k]
                for field in ['peptide','affinity','affinity_percentile','processing_score','presentation_score','presentation_percentile']:
                    assert w[field]==b[field],(mode,k,field)
    counts=[];objectives=[];criterion_summaries=[];crossings=[]
    for mode in ['with_flanks','without_flanks']:
        for scope,fn in [('signal_or_junction',lambda s,e:s<=18),('full_precursor',lambda s,e:True),('unchanged_INS_core',lambda s,e:s>=19),('flank_affected',lambda s,e:s<=22)]:
            for c in ids:
                for model in ['NetMHC','MHCflurry']:
                    rows=raw['NetMHC' if model=='NetMHC' else mode][c]
                    col='percentile_rank' if model=='NetMHC' else 'presentation_percentile'
                    for a in ALLELES:
                        for cut in CUTS:
                            n=sum(float(r[col])<cut for (aa,s,e),r in rows.items() if aa==a and fn(s,e))
                            b=sum(float(r[col])<cut for (aa,s,e),r in raw['NetMHC' if model=='NetMHC' else mode][base].items() if aa==a and fn(s,e))
                            record={'construct':c,'model':model,'context':mode,'scope':scope,'HLA':a,'cutoff_percentile':cut,'baseline_count':b,'variant_count':n,'delta':n-b}
                            counts.append(record)
                            if scope=='signal_or_junction':objectives.append(record)
                    for k,r in rows.items():
                        if c!=candidate or not fn(k[1],k[2]):continue
                        b=raw['NetMHC' if model=='NetMHC' else mode][base][k]
                        for cut in CUTS:
                            if float(r[col])<cut<=float(b[col]):crossings.append({'model':model,'context':mode,'scope':scope,'HLA':k[0],'start':k[1],'end':k[2],'cutoff':cut})
            for c in ids:
                if scope!='signal_or_junction':continue
                rr=[r for r in objectives if r['construct']==c and r['context']==mode]
                assert len(rr)==42
                criterion_summaries.append({'construct':c,'context':mode,'better':sum(r['delta']<0 for r in rr),'unchanged':sum(r['delta']==0 for r in rr),'worse':sum(r['delta']>0 for r in rr)})
    assert not crossings,'A new threshold-crossing pair appeared'
    chosen=[r for r in criterion_summaries if r['construct']==candidate]
    assert [(r['better'],r['unchanged'],r['worse']) for r in chosen]==[(16,26,0),(15,27,0)]
    doubles=[c['id'] for c in constructs if len(c['positions'])==2];assert len(doubles)==8
    feasible=[c for c in doubles if all(r['worse']==0 for r in criterion_summaries if r['construct']==c)]
    assert feasible==[candidate]
    table=[]
    for mode in ['with_flanks','without_flanks']:
        for scope in ['signal_or_junction','full_precursor','unchanged_INS_core','flank_affected']:
            for model in ['NetMHC','MHCflurry']:
                for c in ids:
                    rr=[r for r in counts if r['context']==mode and r['scope']==scope and r['model']==model and r['construct']==c]
                    if rr:table.append({'construct':c,'context':mode,'scope':scope,'model':model,**{f'count_lt_{cut}':sum(r['variant_count'] for r in rr if r['cutoff_percentile']==cut) for cut in CUTS}})
    lookup={(r['construct'],r['context'],r['scope'],r['model']):[r[f'count_lt_{cut}'] for cut in CUTS] for r in table}
    assert lookup[(base,'with_flanks','signal_or_junction','NetMHC')]==[4,8,17]
    assert lookup[(candidate,'with_flanks','signal_or_junction','NetMHC')]==[0,2,8]
    assert lookup[(base,'with_flanks','signal_or_junction','MHCflurry')]==[9,19,27]
    assert lookup[(candidate,'with_flanks','signal_or_junction','MHCflurry')]==[3,6,14]
    for mode in ['with_flanks','without_flanks']:
        assert lookup[(candidate,mode,'unchanged_INS_core','MHCflurry')]==lookup[(base,mode,'unchanged_INS_core','MHCflurry')]
    csv_write(out/'prediction_counts_by_allele.csv',counts);csv_write(out/'objectives42.csv',objectives)
    csv_write(out/'criteria42_summary.csv',criterion_summaries);csv_write(out/'prediction_totals.csv',table)
    expanded=[r for r in counts if r['scope']=='flank_affected']
    expanded_summary=[]
    for mode in ['with_flanks','without_flanks']:
        for c in ids:
            rr=[r for r in expanded if r['construct']==c and r['context']==mode]
            assert len(rr)==42
            expanded_summary.append({'construct':c,'context':mode,'window_HLA_pairs':616,'better':sum(r['delta']<0 for r in rr),'unchanged':sum(r['delta']==0 for r in rr),'worse':sum(r['delta']>0 for r in rr)})
    assert lookup[(base,'with_flanks','flank_affected','NetMHC')]==[4,8,19]
    assert lookup[(candidate,'with_flanks','flank_affected','NetMHC')]==[0,2,10]
    assert lookup[(base,'with_flanks','flank_affected','MHCflurry')]==[9,19,28]
    assert lookup[(candidate,'with_flanks','flank_affected','MHCflurry')]==[3,6,15]
    expanded_chosen=[r for r in expanded_summary if r['construct']==candidate]
    assert [(r['better'],r['unchanged'],r['worse']) for r in expanded_chosen]==[(16,26,0),(15,27,0)]
    csv_write(out/'expanded_objectives42.csv',expanded);csv_write(out/'expanded_criteria42_summary.csv',expanded_summary)
    # One curated table spans every full precursor window, bothMF modes, and the
    # complete NetMHC coverage; origin explicitly distinguishes direct variant
    # predictions from reused exact-sequence canonical scores.
    wide=[]
    for i,c in enumerate(ids,1):
        for k,w in sorted(raw['with_flanks'][c].items()):
            nf=raw['without_flanks'][c][k];el=net['el'][i].get(k);ba=net['ba'][i].get(k)
            r={'construct':c,'HLA':k[0],'start_1based':k[1],'end_1based':k[2],'peptide':w['peptide'],'n_flank_5aa':w['n_flank'],'c_flank_5aa':w['c_flank'],'signal_or_junction':k[1]<=18,'net_score_origin':'direct_variant_context' if k in direct_net_keys[i] else 'exact_peptide_reused_canonical_full_scan','net_EL_percentile':el['percentile_rank'] if el else '', 'net_BA_nM':ba['ic50'] if ba else '', 'net_BA_percentile':ba['percentile_rank'] if ba else '', 'mhc_affinity_nM':w['affinity'],'mhc_affinity_percentile':w['affinity_percentile'],'mhc_presentation_percentile_with_flanks':w['presentation_percentile'],'mhc_presentation_percentile_without_flanks':nf['presentation_percentile'],'mhc_processing_with_flanks':w['processing_score'],'mhc_processing_without_flanks':nf['processing_score']}
            wide.append(r)
    assert len(wide)==34762;csv_write(out/'final_peptide_HLA_predictions.csv',wide)
    reference_validation = reference_comparison(stage_results, out)
    single_validation = single_mhc_analysis(stage_results, out, refs)
    inflammation_projection = inflammation_analysis(out, refs)
    donor=donor_analysis(out,refs)
    infl=read_csv(ROOT/'tables/inflammation_selected_targets.csv');b10=next(r for r in infl if r['peptide']=='HLVEALYLV')
    main=float(b10['mean_ifna_intensity'])/float(b10['mean_basal_intensity']);second=float(b10['independent_median_ifna'])/float(b10['independent_median_basal'])
    assert math.isclose(main,11.252809926850855,rel_tol=1e-12) and math.isclose(second,1.4380281782511148,rel_tol=1e-12)
    novelty=json_read(ROOT/'provenance/novelty/proteome_exact_matches.json')
    assert novelty['protein_sequence_records']==169634
    assert not novelty['matches']['VTAFVNQHL']['exact'] and not novelty['matches']['TAFVNQHL']['IL_equivalent']
    report={'status':'PASS','checksummed_bundle_files':file_count,'runtime':'Python standard library; no network or predictor execution','staged_nonbaseline_variants':{'leader_positions3_to9':133,'position17':19,'selected_doubles':8},'only_tested_double_no_worse_both_contexts':feasible,'reference_comparison':reference_validation,'single_scan_validation':single_validation,'final_constructs':13,'final_prediction_rows':34762,'NetMHC_rows_final_per_output':7462,'MHCflurry_rows_final_per_context':34762,'signal_junction_window_HLA_pairs_per_construct':504,'flank_affected_pairs_per_construct':616,'NetMHC_pairs_reused_after_exact_identity':reused_net_pairs,'L9C_T17C_expanded_criteria':expanded_chosen,'full_precursor_window_HLA_pairs_per_construct':2674,'L9C_T17C_criteria':chosen,'new_pairs_crossing_any_specified_cutoff':0,'donor_projection':donor,'inflammation_projection':inflammation_projection,'B10_18_IFNa_fold_main':main,'B10_18_IFNa_fold_second':second,'sequence_novelty_status':'Cached audited result checked; full excluded proteome not rescanned by this command','source_extraction_status':'Uses attributed target-level workbook excerpts; original workbooks not redistributed','prediction_status':'Cached outputs, not regenerated neural-network predictions','all_mature_products_equal_INS25_110':True}
    (out/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
