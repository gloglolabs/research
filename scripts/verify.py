#!/usr/bin/env python3
"""Verify registered study packages and reproduce their reference CSV tables."""

import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    manifests = sorted((ROOT / "studies").glob("*/*/study.json"))
    if not manifests:
        raise RuntimeError("No registered studies found")
    for manifest in manifests:
        study = json.loads(manifest.read_text())
        config = study["reproduction"]
        base = manifest.parent
        with tempfile.TemporaryDirectory(prefix="gloglo-research-") as scratch:
            output = Path(scratch) / "results"
            subprocess.run(
                [sys.executable, "-I", str(base / config["entrypoint"]),
                 "--output", str(output)],
                cwd=scratch, check=True, capture_output=True, text=True,
            )
            report = json.loads((output / "validation.json").read_text())
            if report.get("status") != "PASS":
                raise RuntimeError(f"{study['id']}: analysis validation failed")
            expected = set(config["expected_csv_files"])
            actual = {path.name for path in output.glob("*.csv")}
            if actual != expected:
                raise RuntimeError(f"{study['id']}: unexpected output inventory: {actual ^ expected}")
            reference = base / config["reference_results"]
            for name in sorted(expected):
                if digest(output / name) != digest(reference / name):
                    raise RuntimeError(f"{study['id']}: reference result differs: {name}")
            print(f"PASS {study['id']}: {len(expected)} identical CSV tables; "
                  f"{report['checksummed_bundle_files']} input checksums verified")


if __name__ == "__main__":
    try:
        main()
    except subprocess.CalledProcessError as error:
        print(error.stdout, file=sys.stderr)
        print(error.stderr, file=sys.stderr)
        raise SystemExit(error.returncode) from error
