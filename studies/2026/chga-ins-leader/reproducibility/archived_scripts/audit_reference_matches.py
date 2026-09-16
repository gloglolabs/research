"""Exact and I/L-equivalent matching against both UniProt human reference files."""
import gzip
import hashlib
import json
import xml.etree.ElementTree as ET
from pathlib import Path

P = Path(__file__).resolve().parent
PEPTIDES = ["VTAFVNQHL", "TAFVNQHL", "VGAFVNQHL", "GAFVNQHL", "VCAFVNQHL", "CAFVNQHL", "VDAFVNQHL", "DAFVNQHL", "AAAFVNQHL", "RSAAVLALL"]
FILES = ["UP000005640_9606.fasta.gz", "UP000005640_9606_additional.fasta.gz"]


def fasta_records(stream):
    header, parts = None, []
    for line in stream:
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(parts)
            header, parts = line[1:].strip(), []
        else:
            parts.append(line.strip())
    if header is not None:
        yield header, "".join(parts)


def main():
    ns = {"m": "http://www.metalinker.org/"}
    release = ET.parse(P / "RELEASE.metalink").getroot()
    files = {e.attrib["name"]: e for e in release.findall("m:files/m:file", ns)}
    records = []
    provenance = []
    for name in FILES:
        data = (P / name).read_bytes()
        expected_md5 = files[name].find("m:verification/m:hash", ns).text
        actual_md5 = hashlib.md5(data).hexdigest()
        assert expected_md5 == actual_md5
        with gzip.open(P / name, "rt") as f:
            batch = list(fasta_records(f))
        records.extend(batch)
        provenance.append({"file": name, "published_md5": expected_md5, "verified": True, "sha256": hashlib.sha256(data).hexdigest(), "sequence_count": len(batch)})
    result = {peptide: {"exact": [], "IL_equivalent": []} for peptide in PEPTIDES}
    for header, sequence in records:
        collapsed = sequence.replace("I", "L")
        for peptide in PEPTIDES:
            for label, target, haystack in [("exact", peptide, sequence), ("IL_equivalent", peptide.replace("I", "L"), collapsed)]:
                offset = haystack.find(target)
                while offset >= 0:
                    result[peptide][label].append({"protein": header, "position_1based": offset + 1, "actual_sequence": sequence[offset:offset+len(peptide)]})
                    offset = haystack.find(target, offset + 1)
    summary = {"UniProt_release": release.find("m:version", ns).text, "reference_proteome": "UP000005640", "protein_sequence_records": len(records), "total_residues": sum(len(s) for _, s in records), "distinct_accessions": len({h.split()[0] for h, _ in records}), "source_files": provenance, "matches": result, "scope": "Exact sequence representation in the specified human reference FASTA release, including its additional sequences. This does not establish T-cell novelty, cross-reactivity, or absence from all noncanonical/variant human peptides."}
    (P / "proteome_exact_matches.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps({"release": summary["UniProt_release"], "records": len(records), "counts": {p: {label: len(v) for label, v in r.items()} for p, r in result.items()}}, indent=2))


if __name__ == "__main__":
    main()
