"""Exploratory single-residue sensitivity at CHGA leader residue17.

Not a functional construct design. All candidates keep CHGA residues1-16,
CHGA18 and the complete native INS25-110 sequence fixed. This explores
whether model-predicted junction binding is sensitive to the leader's -2
residue. Secretion, cleavage, and antigen presentation are not established.
"""
import concurrent.futures
import hashlib
import json
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE = HERE.parent.parent
ENDPOINT = "https://tools-cluster-interface.iedb.org/tools_api/mhci/"
ALLELES = ["HLA-A*02:01", "HLA-A*24:02", "HLA-A*03:01", "HLA-B*40:01", "HLA-B*49:01", "HLA-C*03:04", "HLA-C*07:01"]
LENGTHS = [8, 9, 10, 11]


def main():
    ins = json.loads((BASE / "data/uniprot_P01308.json").read_text())["sequence"]["value"]
    chga = json.loads((BASE / "data/uniprot_P10645.json").read_text())["sequence"]["value"]
    assert chga[16:18] == "TA"
    # 11 downstream residues suffice to cover every8-11mer overlapping residue17.
    contexts = [{"id": "CHGA_T17" + aa, "substitution": "T17" + aa, "sequence": chga[:16] + aa + chga[17:18] + ins[24:35], "variable_position_1based": 17, "junction_after_1based": 18} for aa in "ACDEFGHIKLMNPQRSTVWY"]
    fasta = "\n".join(">" + c["id"] + "\n" + c["sequence"] for c in contexts)
    manifest = {"rationale": __doc__, "alleles": ALLELES, "lengths": LENGTHS, "construct_contexts": contexts, "created_at_utc": datetime.now(timezone.utc).isoformat()}
    (HERE / "variant_manifest.json").write_text(json.dumps(manifest, indent=2))
    (HERE / "variant_contexts.fasta").write_text(fasta + "\n")

    def run(method):
        result_path = HERE / ("variants_" + method + ".tsv")
        if result_path.exists():
            return {"method": method, "status": "cached"}
        parameters = {"method": method, "sequence_text": fasta, "allele": ",".join(a for a in ALLELES for _ in LENGTHS), "length": ",".join(str(n) for _ in ALLELES for n in LENGTHS)}
        payload = urllib.parse.urlencode(parameters).encode()
        (HERE / ("variants_" + method + ".request.json")).write_text(json.dumps({"endpoint": ENDPOINT, "parameters": parameters, "submitted_utc": datetime.now(timezone.utc).isoformat(), "payload_sha256": hashlib.sha256(payload).hexdigest()}, indent=2))
        try:
            with urllib.request.urlopen(urllib.request.Request(ENDPOINT, payload), timeout=180) as response:
                content = response.read()
            result_path.write_bytes(content)
            status = {"method": method, "status": "completed", "bytes": len(content), "sha256": hashlib.sha256(content).hexdigest()}
        except Exception as exc:
            status = {"method": method, "status": "failed", "exception": str(exc)}
        (HERE / ("variants_" + method + ".status.json")).write_text(json.dumps(status, indent=2))
        return status

    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        for result in pool.map(run, ["netmhcpan_ba-4.1", "netmhcpan_el-4.1"]):
            print(json.dumps(result), flush=True)


if __name__ == "__main__":
    main()
