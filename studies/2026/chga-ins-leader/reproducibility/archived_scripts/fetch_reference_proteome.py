"""Fetch the published UniProt reference FASTA files rather than dynamic streaming."""
import concurrent.futures
import gzip
import hashlib
import json
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

P = Path(__file__).resolve().parent
BASE = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/"
FILES = ["UP000005640_9606.fasta.gz", "UP000005640_9606_additional.fasta.gz", "RELEASE.metalink"]


def fetch(name):
    target = P / name
    if target.exists():
        print("Already downloaded", name, flush=True)
        return
    req = urllib.request.Request(BASE + name, headers={"User-Agent": "SequenceAudit/1.0"})
    with urllib.request.urlopen(req, timeout=60) as response:
        headers = dict(response.headers)
        data = response.read()
    if name.endswith(".gz"):
        gzip.decompress(data)  # Verify completion and gzip CRC before naming valid output.
    target.write_bytes(data)
    manifest = {"url": BASE + name, "retrieved_at_utc": datetime.now(timezone.utc).isoformat(), "headers": headers, "sha256": hashlib.sha256(data).hexdigest(), "bytes": len(data)}
    target.with_suffix(target.suffix + ".manifest.json").write_text(json.dumps(manifest, indent=2))
    print("Downloaded", name, len(data), flush=True)


if __name__ == "__main__":
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        list(pool.map(fetch, FILES))
