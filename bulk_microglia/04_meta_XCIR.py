#!/usr/bin/env python3
"""XCIR meta-analysis within regions and across one sample per participant."""

import argparse
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

REGIONS = ["MFG", "STG", "SVZ", "THA", "CC", "CER", "HIP", "OCC", "SN"]
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
    beta = pd.to_numeric(df["tau"], errors="coerce")
    var  = pd.to_numeric(df["var_tau"], errors="coerce")
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
    tau = pd.to_numeric(df["tau"], errors="coerce")
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


def meta_analysis(bb: pd.DataFrame, region: str) -> pd.DataFrame:
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


def select_one_sample(bb_all: pd.DataFrame, metadata: pd.DataFrame,
                      vcf_dir: Path, chromosome: str):
    bb_all = bb_all.merge(
        metadata[["sample_id", "participant_id", "cohort"]],
        left_on="sample", right_on="sample_id", how="left",
    ).drop(columns="sample_id")
    # Select by total XCIR-input variants, not separately for each gene.
    variants = {
        row.sample: count_variants_in_vcfgz(
            vcf_dir / row.cohort / f"{row.sample}.{chromosome}.XCIR_input.vcf.gz"
        )
        for row in bb_all[["sample", "cohort"]].drop_duplicates().itertuples(index=False)
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
    parser.add_argument("--metadata", required=True,
                        help="TSV: sample_id, participant_id, cohort, sex, region")
    parser.add_argument("--results-dir", required=True, type=Path)
    parser.add_argument("--vcf-dir", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--prefix", default="MiGA.ALL")
    parser.add_argument("--chromosome", default="chrX")
    args = parser.parse_args()

    metadata = pd.read_csv(args.metadata, sep="\t", dtype=str)
    metadata = metadata.loc[metadata["sex"] == "Female"].copy()
    args.outdir.mkdir(parents=True, exist_ok=True)

    parts = []
    for region in REGIONS:
        samples = metadata.loc[metadata["region"] == region]
        if samples.empty:
            continue
        regional = []
        for sample in samples.itertuples(index=False):
            path = args.results_dir / sample.cohort / sample.sample_id / "betaBinomXI_results.tsv.gz"
            try:
                bb = pd.read_csv(path, sep="\t")
            except FileNotFoundError:
                print(f"Missing XCIR result: {path}")
                continue
            regional.append(bb)
        bb = pd.concat(regional, ignore_index=True)
        bb["region"] = region
        parts.append(bb)
        result = meta_analysis(bb, region)
        result.to_csv(args.outdir / f"{args.prefix}.{region}.XCIR_gene_meta.tsv.gz",
                      sep="\t", index=False)

    bb_all = pd.concat(parts, ignore_index=True)
    bb_f, selected = select_one_sample(bb_all, metadata, args.vcf_dir, args.chromosome)
    result = meta_analysis(bb_f, "ALLREGIONS")
    result.to_csv(args.outdir / f"{args.prefix}.ALLREGIONS.XCIR_gene_meta.tsv.gz",
                  sep="\t", index=False)
    selected.to_csv(args.outdir / f"{args.prefix}.selected_samples.tsv",
                    sep="\t", index=False)
    print(f"ALLREGIONS: {int(result['escape_call'].sum())} / {len(result)}; {len(selected)} participants")


if __name__ == "__main__":
    main()
