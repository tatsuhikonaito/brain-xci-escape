#!/usr/bin/env python3
"""MiGA per-region and donor-deduplicated XCIR meta-analysis.

Run from the repository root; every path in the input TSVs is relative to that
working directory. All sample / gene eligibility must be established upstream.
"""
import argparse
import gzip
import json
import platform
from pathlib import Path
import numpy as np
import pandas as pd
import scipy
from scipy.stats import norm

FDR_THR = 0.05
RHO_THR = 0.10
MIN_SAMPLES_FOR_META = 3


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    q = np.full_like(p, np.nan, dtype=float)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return q
    p_ok = p[ok]
    n = p_ok.size
    order = np.argsort(p_ok)
    ranks = np.arange(1, n + 1, dtype=float)
    q_ok = p_ok[order] * n / ranks
    q_ok = np.minimum.accumulate(q_ok[::-1])[::-1]
    q_ok = np.clip(q_ok, 0.0, 1.0)
    q_tmp = np.empty_like(p_ok)
    q_tmp[order] = q_ok
    q[ok] = q_tmp
    return q

def ivw_meta_beta(df: pd.DataFrame) -> pd.Series:
    beta = pd.to_numeric(df["tau"], errors="raise")
    var  = pd.to_numeric(df["var_tau"], errors="raise")
    ok = beta.notna() & var.notna() & (var > 0)
    beta = beta[ok]
    var  = var[ok]
    n = int(ok.sum())
    if n == 0:
        return pd.Series({"n_meta": 0, "beta_meta": np.nan, "se_meta": np.nan, "z_meta": np.nan, "p_meta": np.nan, "sum_w": 0.0})
    w = 1.0 / var
    sum_w = float(w.sum())
    beta_meta = float((w * beta).sum() / sum_w)
    se_meta = float(np.sqrt(1.0 / sum_w))
    z_meta = float(beta_meta / se_meta) if se_meta > 0 else np.nan
    p_meta = float(2.0 * norm.sf(abs(z_meta))) if np.isfinite(z_meta) else np.nan
    return pd.Series({"n_meta": n, "beta_meta": beta_meta, "se_meta": se_meta, "z_meta": z_meta, "p_meta": p_meta, "sum_w": sum_w})

def rho_from_tau(df: pd.DataFrame) -> pd.Series:
    tau = pd.to_numeric(df["tau"], errors="raise")
    rho = tau
    ok = rho.notna()
    n = int(ok.sum())
    if n == 0:
        return pd.Series({"n_rho": 0, "rho_mean": np.nan, "rho_median": np.nan, "rho_prop_gt_0.10": np.nan})
    r = rho[ok].astype(float)
    return pd.Series({"n_rho": n, "rho_mean": float(r.mean()), "rho_median": float(r.median()), "rho_prop_gt_0.10": float((r > RHO_THR).mean())})

def count_variants_in_vcfgz(path: str) -> int:
    # Count non-header VCF records (lines not starting with '#')
    n = 0
    with gzip.open(path, "rt") as f:
        for line in f:
            if line and line[0] != "#":
                n += 1
    return n


def summarize(bb: pd.DataFrame, region: str) -> pd.DataFrame:
    # Group selection avoids changing pandas' grouping-column behavior across versions.
    grouped = bb.groupby("GENE", sort=False)[["tau", "var_tau"]]
    meta = grouped.apply(ivw_meta_beta).reset_index()
    rho = grouped.apply(rho_from_tau).reset_index()
    gene_df = meta.merge(rho, on="GENE", how="outer")
    gene_df["pass_min_meta"] = gene_df["n_meta"] >= MIN_SAMPLES_FOR_META
    gene_df["padj_meta"] = np.nan
    keep = gene_df["pass_min_meta"]
    gene_df.loc[keep, "padj_meta"] = bh_fdr(gene_df.loc[keep, "p_meta"].values)
    gene_df["escape_call"] = (
        (gene_df["padj_meta"] < FDR_THR)
        & (gene_df["rho_mean"] > RHO_THR)
        & gene_df["pass_min_meta"]
    )
    gene_df.insert(0, "region", region)
    return gene_df.sort_values(
        ["escape_call", "padj_meta", "rho_mean", "n_meta"],
        ascending=[False, True, False, False],
    ).reset_index(drop=True)


def select_one_sample(bb_all: pd.DataFrame, sample_key: pd.DataFrame):
    # The source notebook chooses a single sample per participant, NOT per gene.
    # Retain input order and the original two-column sort for ties.
    bb_all = bb_all.merge(
        sample_key[["sample_id", "participant_id", "vcf"]],
        left_on="sample", right_on="sample_id", how="left", validate="many_to_one",
    ).drop(columns="sample_id")
    if bb_all[["participant_id", "vcf"]].isna().any().any():
        missing = bb_all.loc[bb_all["participant_id"].isna() | bb_all["vcf"].isna(), "sample"].unique()
        raise ValueError(f"Missing participant / VCF mapping: {missing.tolist()}")
    variants = {
        row.sample: count_variants_in_vcfgz(row.vcf)
        for row in bb_all[["sample", "vcf"]].drop_duplicates().itertuples(index=False)
    }
    bb_all["n_variants_vcf"] = bb_all["sample"].map(variants)
    selected = (
        bb_all[["participant_id", "sample", "n_variants_vcf"]]
        .drop_duplicates()
        .sort_values(["participant_id", "n_variants_vcf"], ascending=[True, False])
        .drop_duplicates("participant_id", keep="first")
    )
    bb_f = bb_all.merge(selected[["participant_id", "sample"]],
                        on=["participant_id", "sample"], how="inner")
    return bb_f, selected


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--region-files", required=True, help="TSV columns: region, path")
    parser.add_argument("--sample-key", required=True, help="TSV columns: sample_id, participant_id, vcf")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--prefix", default="MiGA.ALL")
    args = parser.parse_args()
    regions = pd.read_csv(args.region_files, sep="\t", dtype=str)
    key = pd.read_csv(args.sample_key, sep="\t", dtype=str)
    for label, table, required in [
        ("region-files", regions, ["region", "path"]),
        ("sample-key", key, ["sample_id", "participant_id", "vcf"]),
    ]:
        missing = set(required) - set(table.columns)
        if missing:
            raise ValueError(f"Missing {label} columns: {sorted(missing)}")
        if table.empty or table[required].isna().any().any():
            raise ValueError(f"Empty or incomplete {label} table")
    if regions["region"].duplicated().any() or key["sample_id"].duplicated().any():
        raise ValueError("region-files must have unique regions; sample-key must have unique sample IDs")
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    parts = []
    for entry in regions.itertuples(index=False):
        bb = pd.read_csv(entry.path, sep="\t", dtype={"GENE": str, "sample": str})
        required = ["GENE", "sample", "tau", "var_tau"]
        missing = set(required) - set(bb.columns)
        if missing:
            raise ValueError(f"Missing XCIR columns in {entry.path}: {sorted(missing)}")
        if bb.empty or bb[["GENE", "sample"]].isna().any().any():
            raise ValueError(f"Empty result or missing gene / sample IDs in {entry.path}")
        for col in ["tau", "var_tau"]:
            bb[col] = pd.to_numeric(bb[col], errors="raise")
        bb["region"] = entry.region
        parts.append(bb)
        result = summarize(bb, entry.region)
        result.to_csv(outdir / f"{args.prefix}.{entry.region}.XCIR_gene_meta.tsv.gz", sep="\t", index=False)
        print(f"{entry.region}: {int(result['escape_call'].sum())} / {len(result)}")
    bb_all = pd.concat(parts, ignore_index=True)
    bb_f, selected = select_one_sample(bb_all, key)
    result = summarize(bb_f, "ALLREGIONS")
    result.to_csv(outdir / f"{args.prefix}.ALLREGIONS.XCIR_gene_meta.tsv.gz", sep="\t", index=False)
    selected.to_csv(outdir / f"{args.prefix}.selected_samples.tsv", sep="\t", index=False)
    run_info = {"python": platform.python_version(), "numpy": np.__version__,
                "pandas": pd.__version__, "scipy": scipy.__version__,
                "regions": regions["region"].tolist(), "FDR_THR": FDR_THR,
                "RHO_THR": RHO_THR, "MIN_SAMPLES_FOR_META": MIN_SAMPLES_FOR_META,
                "selected_participants": len(selected), **vars(args)}
    (outdir / f"{args.prefix}.run_info.json").write_text(json.dumps(run_info, indent=2) + "\n")
    print(f"ALLREGIONS: {int(result['escape_call'].sum())} / {len(result)}; {len(selected)} participants")


if __name__ == "__main__":
    main()
