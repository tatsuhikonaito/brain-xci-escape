#!/usr/bin/env python3
"""Mechanical collection of an explicitly supplied analysis sample list.

This is a public input adapter, not a reconstruction of the missing historical
sample-selection / per-region assembly script. No sex, skewness, gene annotation
or XCI-state filter is inferred here. All paths refer to the working directory.
"""
import argparse
from pathlib import Path
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True,
                        help="TSV: sample_id, region, result_file; analysis-eligible samples only")
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()
    manifest = pd.read_csv(args.manifest, sep="\t", dtype=str)
    required = ["sample_id", "region", "result_file"]
    missing = set(required) - set(manifest.columns)
    if missing:
        raise ValueError(f"Missing manifest columns: {sorted(missing)}")
    if manifest.empty or manifest[required].isna().any().any():
        raise ValueError("The manifest must contain non-missing sample IDs, regions and paths")
    if manifest.duplicated(["sample_id", "region"]).any():
        raise ValueError("Duplicate sample_id / region entries in manifest")
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    region_files = []
    for region, entries in manifest.groupby("region", sort=False):
        parts = []
        for row in entries.itertuples(index=False):
            bb = pd.read_csv(row.result_file, sep="\t", dtype={"sample": str, "GENE": str})
            if not {"sample", "GENE", "tau", "var_tau"}.issubset(bb.columns):
                raise ValueError(f"Required XCIR columns missing in {row.result_file}")
            if bb.empty:
                raise ValueError(f"Empty XCIR result for manifest sample {row.sample_id}")
            if bb["sample"].isna().any() or not bb["sample"].eq(row.sample_id).all():
                raise ValueError(f"XCIR sample labels do not match {row.sample_id}; supply the correct ID mapping upstream")
            parts.append(bb)
        output = outdir / f"{region}.betaBinomXI_results.tsv.gz"
        pd.concat(parts, ignore_index=True).to_csv(output, sep="\t", index=False)
        region_files.append({"region": region, "path": str(output)})
    pd.DataFrame(region_files).to_csv(outdir / "region_files.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
