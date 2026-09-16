"""Explicit double candidates, including overlap of the two edited positions."""
import json
from scan_leader import HERE, ALLELES, LENGTHS, utc

root = json.loads((HERE / "scan_manifest.json").read_text())
baseline = root["constructs"][0]["full_sequence"]
designs = [[], [(9,"C")], [(9,"G")], [(17,"C")], [(17,"G")]]
designs += [[(9,a), (17,b)] for a in "CG" for b in "CGDP"]
out = HERE / "double_scan"
out.mkdir(parents=True, exist_ok=True)
variants = []
for idx, changes in enumerate(designs,1):
    sequence = list(baseline)
    names = []
    for position, aa in changes:
        names.append(f"{baseline[position-1]}{position}{aa}")
        sequence[position-1] = aa
    sequence = "".join(sequence)
    assert sequence[18:] == baseline[18:]
    assert sum(a != b for a,b in zip(sequence,baseline)) == len(changes)
    name = "_".join(names) if names else "baseline"
    variants.append({"seq_num":idx, "id":"CHGA_"+name, "substitution":name, "positions":[p for p,a in changes],
                     "leader_sequence":sequence[:18], "full_sequence":sequence, "sequence":sequence[:29],
                     "junction_after_1based":18})
manifest = {"scope":__doc__, "selection":"L9C/G are the only leader3-9 single substitutions without leader-only NetMHCpan EL<0.5 hits; combine with priorT17C/G/D/P frontier, not a validated functional design.",
            "alleles":ALLELES,"lengths":LENGTHS,"constructs":variants,"created_utc":utc(),"expected_rows_per_method":len(variants)*574}
(out / "scan_manifest.json").write_text(json.dumps(manifest,indent=2))
(out / "variant_contexts.fasta").write_text("\n".join(">"+v["id"]+"\n"+v["sequence"] for v in variants)+"\n")
(out / "full_constructs.fasta").write_text("\n".join(">"+v["id"]+"\n"+v["full_sequence"] for v in variants)+"\n")
print(f"{len(variants)}constructs ready at{out}")
