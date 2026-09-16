"""Smaller requests after an uninformative HTTP500 on the initial large BA job."""
import argparse
import csv
import hashlib
import io
import json
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from scan_leader import HERE, ENDPOINT, utc


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("method", choices=["netmhcpan_ba-4.1", "netmhcpan_el-4.1"])
    parser.add_argument("--directory", type=Path, default=HERE)
    args = parser.parse_args()
    directory = args.directory
    manifest = json.loads((directory / "scan_manifest.json").read_text())
    out = directory / "batches" / args.method
    out.mkdir(parents=True, exist_ok=True)
    joined = []
    for offset in range(0, len(manifest["constructs"]), 20):
        batch = manifest["constructs"][offset:offset+20]
        name = f"batch_{offset//20+1:02d}"
        fasta = "\n".join(">" + v["id"] + "\n" + v["sequence"] for v in batch) + "\n"
        params = {"method": args.method, "sequence_text": fasta,
                  "allele": ",".join(a for a in manifest["alleles"] for _ in manifest["lengths"]),
                  "length": ",".join(str(k) for _ in manifest["alleles"] for k in manifest["lengths"])}
        payload = urllib.parse.urlencode(params).encode()
        raw = out / (name + ".tsv")
        if not raw.exists():
            (out / (name + ".request.json")).write_text(json.dumps({"endpoint": ENDPOINT, "parameters": params,
                "global_seq_num_offset": offset, "local_to_global_seq_num": {i+1:v["seq_num"] for i,v in enumerate(batch)},
                "submitted_utc": utc(), "payload_sha256": hashlib.sha256(payload).hexdigest()}, indent=2))
            try:
                with urllib.request.urlopen(urllib.request.Request(ENDPOINT, payload), timeout=600) as response:
                    body = response.read()
                    status_code = response.status
                raw.write_bytes(body)
                (out / (name + ".status.json")).write_text(json.dumps({"status": "completed", "http_status": status_code,
                    "bytes": len(body), "sha256": hashlib.sha256(body).hexdigest(), "completed_utc": utc()}, indent=2))
            except urllib.error.HTTPError as error:
                (out / (name + ".error.txt")).write_bytes(error.read())
                (out / (name + ".status.json")).write_text(json.dumps({"status":"HTTP error", "code":error.code, "reason":error.reason}, indent=2))
                raise
        rows = list(csv.DictReader(io.StringIO(raw.read_text()), delimiter="\t"))
        assert len(rows) == len(batch)*7*82, (name, len(rows))
        for row in rows:
            local = int(row["seq_num"])
            assert 1 <= local <= len(batch)
            row["seq_num"] = str(batch[local-1]["seq_num"])
            joined.append(row)
        print(f"{args.method} {name}: validated{len(rows)}rows", flush=True)
    target = directory / (args.method + ".tsv")
    with target.open("w") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(joined[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(joined)
    assert len(joined) == manifest["expected_rows_per_method"]
    (directory / (args.method + ".batched_status.json")).write_text(json.dumps({"status":"completed", "rows":len(joined),
        "sha256":hashlib.sha256(target.read_bytes()).hexdigest(), "completed_utc":utc()}, indent=2))


if __name__ == "__main__":
    main()
