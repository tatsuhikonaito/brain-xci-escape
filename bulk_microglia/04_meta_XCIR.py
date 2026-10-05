#!/usr/bin/env python3

import argparse
import gzip
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import norm

FDR_THR = 0.05
RHO_THR = 0.10
MIN_SAMPLES_FOR_META = 3


def bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    q = np.full_like(p, np.nan, dtype=float)
    ok = np.isfinite(p)
    p_ok = p[ok]
    n = p_ok.size
    if n == 0:
        return q
    order = np.argsort(p_ok)
    q_ok = p_ok[order] * n / np.arange(1, n + 1)
    q_ok = np.minimum.accumulate(q_ok[::-1])[::-1]
    q_ok = np.clip(q_ok, 0, 1)
    tmp = np.empty_like(p_ok)
    tmp[order] = q_ok
    q[ok] = tmp
    return q


def ivw_meta(df):
    beta = pd.to_numeric(df["tau"], errors="coerce")
    var = pd.to_numeric(df["var_tau"], errors="coerce")
    ok = beta.notna() & var.notna() & (var > 0)
    beta = beta[ok]
    var = var[ok]
    n = int(ok.sum())
    if n == 0:
        return pd.Series({"n_meta": 0, "beta_meta": np.nan, "se_meta": np.nan,
                          "z_meta": np.nan, "p_meta": np.nan, "sum_w": 0.0})
    w = 1 / var
    sum_w = float(w.sum())
    beta_meta = float((w * beta).sum() / sum_w)
    se_meta = float(np.sqrt(1 / sum_w))
    z_meta = beta_meta / se_meta
    p_meta = float(2 * norm.sf(abs(z_meta)))
    return pd.Series({"n_meta": n, "beta_meta": beta_meta, "se_meta": se_meta,
                      "z_meta": z_meta, "p_meta": p_meta, "sum_w": sum_w})


def rho_summary(df):
    tau = pd.to_numeric(df["tau"], errors="coerce").dropna()
    if len(tau) == 0:
        return pd.Series({"n_rho": 0, "rho_mean": np.nan, "rho_median": np.nan,
                          "rho_prop_gt_0.10": np.nan})
    return pd.Series({
        "n_rho": len(tau),
        "rho_mean": float(tau.mean()),
        "rho_median": float(tau.median()),
        "rho_prop_gt_0.10": float((tau > RHO_THR).mean()),
    })


def summarize(bb, label):
    grouped = bb.groupby("GENE", sort=False)[["tau", "var_tau"]]
    meta = grouped.apply(ivw_meta).reset_index()
    rho = grouped.apply(rho_summary).reset_index()
    out = meta.merge(rho, on="GENE", how="outer")
    out["pass_min_meta"] = out["n_meta"] >= MIN_SAMPLES_FOR_META
    out["padj_meta"] = np.nan
    keep = out["pass_min_meta"]
    out.loc[keep, "padj_meta"] = bh_fdr(out.loc[keep, "p_meta"].values)
    out["escape_call"] = (
        (out["padj_meta"] < FDR_THR)
        & (out["rho_mean"] > RHO_THR)
        & out["pass_min_meta"]
    )
    out.insert(0, "region", label)
    return out.sort_values(
        ["escape_call", "padj_meta", "rho_mean", "n_meta"],
        ascending=[False, True, False, False],
    ).reset_index(drop=True)


def count_variants(path):
    n = 0
    with gzip.open(path, "rt") as f:
        for line in f:
            if line and line[0] != "#":
                n += 1
    return n


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True,
                        help="TSV columns: sample_id, participant_id, region, result_file, vcf")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--prefix", default="bulk_microglia")
    args = parser.parse_args()

    manifest = pd.read_csv(args.manifest, sep="\t", dtype=str)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    region_tables = []
    all_tables = []

    for region, m in manifest.groupby("region", sort=False):
        parts = []
        for row in m.itertuples(index=False):
            bb = pd.read_csv(row.result_file, sep="\t", dtype={"GENE": str, "sample": str})
            parts.append(bb)
        bb_region = pd.concat(parts, ignore_index=True)
        all_tables.append(bb_region.assign(region=region))
        result = summarize(bb_region, region)
        result.to_csv(outdir / f"{args.prefix}.{region}.XCIR_gene_meta.tsv.gz",
                      sep="\t", index=False)

    bb_all = pd.concat(all_tables, ignore_index=True)
    sample_key = manifest[["sample_id", "participant_id", "vcf"]].drop_duplicates()
    bb_all = bb_all.merge(sample_key, left_on="sample", right_on="sample_id", how="left")

    variant_counts = {
        row.sample_id: count_variants(row.vcf)
        for row in sample_key.itertuples(index=False)
    }
    bb_all["n_variants_vcf"] = bb_all["sample"].map(variant_counts)

    selected = (
        bb_all[["participant_id", "sample", "n_variants_vcf"]]
        .drop_duplicates()
        .sort_values(["participant_id", "n_variants_vcf"], ascending=[True, False])
        .drop_duplicates("participant_id", keep="first")
    )

    bb_all = bb_all.merge(selected[["participant_id", "sample"]],
                          on=["participant_id", "sample"], how="inner")
    result = summarize(bb_all, "ALLREGIONS")
    result.to_csv(outdir / f"{args.prefix}.ALLREGIONS.XCIR_gene_meta.tsv.gz",
                  sep="\t", index=False)
    selected.to_csv(outdir / f"{args.prefix}.selected_samples.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
