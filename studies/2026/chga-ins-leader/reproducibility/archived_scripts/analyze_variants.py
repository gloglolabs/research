"""Summarize exhaustive CHGA T17 model sensitivity, with paired controls."""
import csv
import json
from pathlib import Path

P = Path(__file__).resolve().parent


def read(name):
    return list(csv.DictReader((P / name).open(), delimiter="\t"))


def key(row):
    return tuple(row[k] for k in ["seq_num", "allele", "start", "end", "peptide"])


def main():
    manifest = json.loads((P / "variant_manifest.json").read_text())
    contexts = manifest["construct_contexts"]
    ba, el = read("variants_netmhcpan_ba-4.1.tsv"), read("variants_netmhcpan_el-4.1.tsv")
    expected = sum(sum(len(c["sequence"]) - n + 1 for n in manifest["lengths"]) for c in contexts) * len(manifest["alleles"])
    assert len(ba) == len(el) == expected
    assert len({key(r) for r in ba}) == expected
    ba_index = {key(r): r for r in ba}
    output, comparisons = [], []
    for i, c in enumerate(contexts, 1):
        rows = [r for r in el if int(r["seq_num"]) == i]
        for r in rows:
            assert c["sequence"][int(r["start"])-1:int(r["end"])] == r["peptide"]
            assert len(r["peptide"]) == int(r["length"])
            assert key(r) in ba_index
        affected = [r for r in rows if int(r["start"]) <= 17 <= int(r["end"])]
        record = {"substitution": c["substitution"], "affected_window_HLA_pairs": len(affected)}
        for cutoff in [.5, 1, 2, 5]:
            record[f"EL_rank_lt_{cutoff}"] = sum(float(r["percentile_rank"]) < cutoff for r in affected)
        record["min_EL_percentile"] = min(float(r["percentile_rank"]) for r in affected)
        record["BA_nM_lt_500"] = sum(float(ba_index[key(r)]["ic50"]) < 500 for r in affected)
        for allele in manifest["alleles"]:
            ar = [r for r in affected if r["allele"] == allele]
            record[allele + "_EL_lt2"] = sum(float(r["percentile_rank"]) < 2 for r in ar)
        output.append(record)
        if c["substitution"] in ["T17T", "T17G", "T17C", "T17S"]:
            for r in affected:
                pair = ba_index[key(r)]
                if float(r["percentile_rank"]) < 2 or float(pair["ic50"]) < 500:
                    comparisons.append({"substitution": c["substitution"], "peptide": r["peptide"], "allele": r["allele"], "EL_percentile_rank": float(r["percentile_rank"]), "BA_percentile_rank": float(pair["percentile_rank"]), "BA_nM": float(pair["ic50"])})
    for name, data in [("variant_summary.csv", output), ("variant_candidate_pairs.csv", comparisons)]:
        with (P / name).open("w") as f:
            writer = csv.DictWriter(f, fieldnames=list(data[0]))
            writer.writeheader()
            writer.writerows(data)
    # The T17T context is an independently submitted control for the original scan.
    reference = list(csv.DictReader((P.parent / "predictions/netmhcpan_el-4.1.tsv").open(), delimiter="\t"))
    control = {tuple(r[k] for k in ["allele", "start", "end", "peptide"]): r for r in el if int(r["seq_num"]) == 17}
    matched = 0
    for r in reference:
        k = tuple(r[x] for x in ["allele", "start", "end", "peptide"])
        if r["seq_num"] == "2" and k in control:
            assert r["score"] == control[k]["score"]
            assert r["percentile_rank"] == control[k]["percentile_rank"]
            matched += 1
    assert matched == expected // len(contexts)
    validation = {"variants_including_baseline": len(contexts), "values_per_model": expected, "combined_model_values": 2 * expected, "baseline_control_matches": matched, "threshold_rule": "Strict less-than on API rounded percentile output; overlapping windows are not independent epitopes.", "scope": "Only windows containing residue17 are counted in variant comparisons. Other CHGA-derived candidate ligands remain unchanged.", "functional_status": "No functional or immunological validation. SignalP results analyzed separately."}
    (P / "variant_validation.json").write_text(json.dumps(validation, indent=2))
    print(json.dumps(validation, indent=2))


if __name__ == "__main__":
    main()
