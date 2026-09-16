"""Exhaustive single-substitution HLA binding scan of canonical CHGA leader3-9.

Canonical baseline plus 133 single substitutions. Whole native INS25-110 remains
unchanged in every inferred construct. 29-residue request contexts cover every
8-11mer in the leader or crossing the leader/INS boundary. No functional claim.
"""
from pathlib import Path
import concurrent.futures
import datetime
import hashlib
import json
import urllib.error
import urllib.parse
import urllib.request

HERE = Path(__file__).resolve().parent
BASE = HERE.parent.parent
ENDPOINT = "https://tools-cluster-interface.iedb.org/tools_api/mhci/"
ALLELES = ["HLA-A*02:01", "HLA-A*24:02", "HLA-A*03:01", "HLA-B*40:01", "HLA-B*49:01", "HLA-C*03:04", "HLA-C*07:01"]
LENGTHS = [8, 9, 10, 11]
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"


def utc():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def main():
    HERE.mkdir(parents=True, exist_ok=True)
    ins = json.loads((BASE / "data/uniprot_P01308.json").read_text())["sequence"]["value"]
    chga = json.loads((BASE / "data/uniprot_P10645.json").read_text())["sequence"]["value"]
    leader = chga[:18]
    full = leader + ins[24:]
    assert leader == "MRSAAVLALLLCAGQVTA"
    reference = json.loads((HERE.parent / "predictions/scan_manifest.json").read_text())
    assert full == reference["constructs"][1]["sequence"]
    variants = [{"id": "CANONICAL_CHGA18_INS25_BASELINE", "substitution": "baseline", "position": None,
                 "native_aa": None, "replacement_aa": None, "leader_sequence": leader,
                 "full_sequence": full, "sequence": full[:29], "junction_after_1based": 18}]
    for pos in range(3, 10):
        for aa in AMINO_ACIDS:
            native = leader[pos-1]
            if aa == native:
                continue
            mutant = leader[:pos-1] + aa + leader[pos:]
            newfull = mutant + ins[24:]
            assert sum(a != b for a, b in zip(full, newfull)) == 1
            assert newfull[18:] == ins[24:]
            variants.append({"id": f"CHGA_{native}{pos}{aa}", "substitution": f"{native}{pos}{aa}",
                             "position": pos, "native_aa": native, "replacement_aa": aa,
                             "leader_sequence": mutant, "full_sequence": newfull, "sequence": newfull[:29],
                             "junction_after_1based": 18})
    assert len(variants) == 134
    assert len({v["sequence"] for v in variants}) == 134
    for i, variant in enumerate(variants, 1):
        variant["seq_num"] = i
    manifest = {"scope": __doc__, "source_manifest": "../predictions/scan_manifest.json",
                "alleles": ALLELES, "lengths": LENGTHS, "positions": list(range(3, 10)),
                "constructs": variants, "expected_rows_per_method": 134 * 7 * sum(29-k+1 for k in LENGTHS),
                "created_utc": utc()}
    (HERE / "scan_manifest.json").write_text(json.dumps(manifest, indent=2))
    fasta = "\n".join(">" + v["id"] + "\n" + v["sequence"] for v in variants) + "\n"
    (HERE / "variant_contexts.fasta").write_text(fasta)
    (HERE / "full_constructs.fasta").write_text("\n".join(">" + v["id"] + "\n" + v["full_sequence"] for v in variants) + "\n")
    print(f"Manifest ready: {len(variants)} constructs, {manifest['expected_rows_per_method']} expected rows per method", flush=True)

    def run(method):
        path = HERE / (method + ".tsv")
        if path.exists():
            return {"method": method, "status": "cached", "bytes": path.stat().st_size}
        params = {"method": method, "sequence_text": fasta,
                  "allele": ",".join(a for a in ALLELES for _ in LENGTHS),
                  "length": ",".join(str(k) for _ in ALLELES for k in LENGTHS)}
        payload = urllib.parse.urlencode(params).encode()
        (HERE / (method + ".request.json")).write_text(json.dumps({"endpoint": ENDPOINT, "parameters": params,
              "submitted_utc": utc(), "payload_sha256": hashlib.sha256(payload).hexdigest()}, indent=2))
        try:
            with urllib.request.urlopen(urllib.request.Request(ENDPOINT, payload), timeout=600) as response:
                body = response.read()
                http_status = response.status
            path.write_bytes(body)
            status = {"method": method, "status": "completed", "http_status": http_status,
                      "bytes": len(body), "sha256": hashlib.sha256(body).hexdigest(), "completed_utc": utc()}
        except urllib.error.HTTPError as error:
            (HERE / (method + ".error.txt")).write_bytes(error.read())
            status = {"method": method, "status": "HTTP error", "code": error.code, "reason": error.reason}
        except (urllib.error.URLError, TimeoutError) as error:
            status = {"method": method, "status": "network error", "reason": str(error)}
        (HERE / (method + ".status.json")).write_text(json.dumps(status, indent=2))
        return status

    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        for result in pool.map(run, ["netmhcpan_ba-4.1", "netmhcpan_el-4.1"]):
            print(json.dumps(result), flush=True)


if __name__ == "__main__":
    main()
