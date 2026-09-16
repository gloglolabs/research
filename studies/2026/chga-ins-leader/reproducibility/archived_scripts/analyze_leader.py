"""Validate paired prediction grids, count every window, and solve hit-set bounds."""
import argparse
import csv
import hashlib
import itertools
import json
from collections import defaultdict
from pathlib import Path
from scan_leader import HERE


def read(path):
    return list(csv.DictReader(path.open(), delimiter="\t"))


def write(path, rows):
    if not rows:
        path.write_text("")
        return
    with path.open("w") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def key(row):
    return row["allele"], int(row["seq_num"]), int(row["start"]), int(row["end"])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--directory", type=Path, default=HERE)
    args = parser.parse_args()
    directory = args.directory
    manifest = json.loads((directory / "scan_manifest.json").read_text())
    variants = {v["seq_num"]:v for v in manifest["constructs"]}
    methods = {m:{key(r):r for r in read(directory / f"netmhcpan_{m}-4.1.tsv")} for m in ("ba","el")}
    expected = {(a,i,s,s+k-1) for a in manifest["alleles"] for i in variants
                for k in manifest["lengths"] for s in range(1,len(variants[i]["sequence"])-k+2)}
    assert len(expected) == manifest["expected_rows_per_method"]
    for method in methods:
        assert len(read(directory / f"netmhcpan_{method}-4.1.tsv")) == len(expected)
        assert set(methods[method]) == expected
    prior = {m:{key(r):r for r in read(HERE.parent / "predictions" / f"netmhcpan_{m}-4.1.tsv")} for m in methods}
    baseline_checks = 0
    all_rows = []
    for a,i,s,e in sorted(expected):
        ba,el = methods["ba"][(a,i,s,e)],methods["el"][(a,i,s,e)]
        peptide = variants[i]["sequence"][s-1:e]
        assert ba["peptide"] == el["peptide"] == peptide
        assert int(ba["length"]) == int(el["length"]) == e-s+1
        base_el = methods["el"][(a,1,s,e)]
        base_ba = methods["ba"][(a,1,s,e)]
        if i == 1:
            for m,r in (("ba",ba),("el",el)):
                original = prior[m][(a,2,s,e)]
                assert r["peptide"] == original["peptide"]
                assert float(r["percentile_rank"]) == float(original["percentile_rank"])
                if m == "ba":
                    assert float(r["ic50"]) == float(original["ic50"])
                baseline_checks += 1
        if peptide == base_el["peptide"]:
            assert float(el["percentile_rank"]) == float(base_el["percentile_rank"])
            assert float(ba["ic50"]) == float(base_ba["ic50"])
        region = "leader" if e <= 18 else "junction" if s <= 18 else "downstream"
        all_rows.append({"variant":variants[i]["substitution"],"seq_num":i,"allele":a,"start":s,"end":e,
            "length":e-s+1,"peptide":peptide,"region":region,"baseline_peptide":base_el["peptide"],
            "peptide_changed":peptide != base_el["peptide"],"EL_rank":float(el["percentile_rank"]),
            "BA_rank":float(ba["percentile_rank"]),"BA_nM":float(ba["ic50"]),
            "baseline_EL_rank":float(base_el["percentile_rank"]),"baseline_BA_nM":float(base_ba["ic50"])})
    grouped = defaultdict(list)
    for row in all_rows:
        grouped[row["seq_num"]].append(row)
    metrics = []
    allele_metrics = []
    for i,v in variants.items():
        metric = {"variant":v["substitution"],"seq_num":i,"leader_sequence":v["leader_sequence"]}
        for region in ("leader","junction","leader_and_junction"):
            subset = [r for r in grouped[i] if r["region"] == region or region == "leader_and_junction" and r["region"] != "downstream"]
            metric[f"{region}_windows"] = len(subset)
            for threshold in (.5,1,2):
                metric[f"{region}_EL_lt{threshold}"] = sum(r["EL_rank"]<threshold for r in subset)
                metric[f"{region}_emergent_EL_lt{threshold}"] = sum(r["EL_rank"]<threshold<=r["baseline_EL_rank"] for r in subset)
                metric[f"{region}_removed_EL_lt{threshold}"] = sum(r["baseline_EL_rank"]<threshold<=r["EL_rank"] for r in subset)
            for threshold in (50,500):
                metric[f"{region}_BA_lt{threshold}nM"] = sum(r["BA_nM"]<threshold for r in subset)
            for allele in manifest["alleles"]:
                ar = [r for r in subset if r["allele"] == allele]
                allele_metrics.append({"variant":v["substitution"],"region":region,"allele":allele,
                    "EL_lt0.5":sum(r["EL_rank"]<.5 for r in ar),"EL_lt1":sum(r["EL_rank"]<1 for r in ar),
                    "EL_lt2":sum(r["EL_rank"]<2 for r in ar),"BA_lt500nM":sum(r["BA_nM"]<500 for r in ar)})
        metrics.append(metric)
    objectives = [f"leader_EL_lt{t}" for t in (.5,1,2)] + ["leader_BA_lt500nM"]
    objectives_all = [f"leader_and_junction_EL_lt{t}" for t in (.5,1,2)] + ["leader_and_junction_BA_lt500nM"]
    for row in metrics:
        for name,columns in (("leader_count_pareto",objectives),("all_count_pareto",objectives_all)):
            row[name] = not any(all(other[c]<=row[c] for c in columns) and any(other[c]<row[c] for c in columns) for other in metrics)
    write(directory / "variant_metrics.csv",metrics)
    write(directory / "allele_metrics.csv",allele_metrics)
    write(directory / "candidate_windows.csv",[r for r in all_rows if r["region"] != "downstream" and (r["EL_rank"]<2 or r["BA_nM"]<500 or r["baseline_EL_rank"]<.5)])
    write(directory / "all_scores.csv",all_rows)
    bounds = []
    for cutoff in (.5,1,2):
        baseline_rows = [r for r in grouped[1] if r["region"] != "downstream" and r["EL_rank"]<cutoff]
        windows = [set(range(r["start"],min(r["end"],18)+1)) for r in baseline_rows]
        minimum_sets = []
        for k in range(1,19):
            minimum_sets = [c for c in itertools.combinations(range(1,19),k) if all(set(c)&w for w in windows)]
            if minimum_sets:
                break
        bounds.append({"EL_strict_cutoff":cutoff,"original_candidate_pairs":len(baseline_rows),
            "mutable_positions":"canonical CHGA leader1-18 only", "minimum_edited_positions":k,
            "all_minimum_position_sets":minimum_sets,
            "interpretation":"Necessary sequence-hit condition only. Does not assert a replacement removes binding or avoids creating other binders."})
    (directory / "minimum_hitsets.json").write_text(json.dumps(bounds,indent=2))
    validation = {"rows_per_method":len(expected),"constructs":len(variants),"complete_grid":True,
        "all_sequences_match_coordinates":True,"unchanged_peptides_have_unchanged_scores":True,
        "baseline_comparisons_to_original_full_construct":baseline_checks,
        "baseline_matches_original_full_construct":True,
        "input_sha256":{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in directory.glob("netmhcpan_*.tsv")},
        "leader_pareto_objectives":objectives,"all_pareto_objectives":objectives_all}
    (directory / "validation.json").write_text(json.dumps(validation,indent=2))
    for row in metrics:
        if row["variant"] == "baseline" or row["leader_count_pareto"] or row["all_count_pareto"] or len(variants)<20:
            print(row["variant"],"leader",[row[c] for c in objectives],"all",[row[c] for c in objectives_all],flush=True)


if __name__ == "__main__":
    main()
