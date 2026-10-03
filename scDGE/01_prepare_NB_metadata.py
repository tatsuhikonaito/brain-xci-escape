#!/usr/bin/env python3
"""Prepare per-cell covariates for negative-binomial regression."""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import harmonypy as hm
import scipy.sparse as sp


def main():
    # Input and output paths (run from the repository root).
    adata_path = "/path/to/raw_counts.h5ad"
    clinical_path = "/path/to/clinical.csv"
    out_path = "/path/to/scDGE/metadata.txt.gz"
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    # Read raw counts.
    if not os.path.exists(adata_path):
        raise FileNotFoundError(f"h5ad file not found: {adata_path}")

    adata = sc.read_h5ad(adata_path)

    # Total UMI from raw counts (for NB covariate)
    X = adata.X
    if sp.issparse(X):
        nUMI = np.asarray(X.sum(axis=1)).ravel()
    else:
        nUMI = X.sum(axis=1)

    # log10(nUMI + 1)
    adata.obs["nUMI"] = np.log10(nUMI + 1.0)

    # Log-normalization for PCA / Harmony
    # library-size normalization
    sc.pp.normalize_total(adata, target_sum=1e4)
    # log-transform
    sc.pp.log1p(adata)

    # PCA & Harmony
    sc.pp.pca(adata, n_comps=30)

    # Run Harmony using individual_ID as batch covariate
    if "individual_ID" not in adata.obs.columns:
        raise KeyError("adata.obs does not contain 'individual_ID' column.")

    harmony = hm.run_harmony(
        adata.obsm["X_pca"],
        adata.obs,
        ["individual_ID"],
    )
    adata.obsm["X_pca_harmony"] = harmony.Z_corr.T

    # Build per-cell metadata (projid / individualID / nUMI / Annotation / PCs)
    required_cols = ["projid", "individual_ID", "nUMI", "major_cell_type"]
    for col in required_cols:
        if col not in adata.obs.columns:
            raise KeyError(f"adata.obs does not contain required column: {col}")

    metadata = adata.obs[required_cols].copy()
    metadata = metadata.rename(
        columns={
            "individual_ID": "individualID",
            "major_cell_type": "Annotation",
        }
    )

    # original PCs
    pca_df = pd.DataFrame(
        adata.obsm["X_pca"][:, :10],
        index=adata.obs.index,
        columns=[f"oriPC_{i + 1}" for i in range(10)],
    )

    # Harmony PCs
    hm_pca_df = pd.DataFrame(
        adata.obsm["X_pca_harmony"][:, :10],
        index=adata.obs.index,
        columns=[f"hmPC_{i + 1}" for i in range(10)],
    )

    # combine
    metadata = pd.concat([metadata, pca_df, hm_pca_df], axis=1)

    # index = barcode
    metadata.index.name = "barcode"
    metadata = metadata.reset_index()

    # Merge with clinical metadata
    if not os.path.exists(clinical_path):
        raise FileNotFoundError(f"Clinical file not found: {clinical_path}")

    clinical = pd.read_csv(clinical_path)

    # APOE4 carrier flag
    clinical["apoe4"] = np.where(
        clinical["apoe_genotype"].isna(),
        np.nan,
        np.where(
            clinical["apoe_genotype"].astype(str).str.contains("4"),
            1,
            0,
        ),
    )

    # Avoid duplicate individualID column before merge
    if "individualID" in metadata.columns:
        metadata = metadata.drop(columns=["individualID"])

    merge_cols = [
        "projid",
        "individualID",
        "age_death",
        "cogdx",
        "pmi",
        "msex",
        "Study",
        "apoe4",
    ]
    missing_cols = [c for c in merge_cols if c not in clinical.columns]
    if missing_cols:
        raise KeyError(
            "Clinical file is missing required columns: " + ", ".join(missing_cols)
        )

    metadata = (
        metadata.merge(
            clinical[merge_cols],
            on="projid",
            how="left",
        )
        .assign(
            age=lambda df: df["age_death"].replace("90+", "90").astype(float),
            cogdx=lambda df: df["cogdx"].astype(str),
            individualID=lambda df: df["individualID"].astype(str),
            # Convert msex: 0 -> 1 (female), 1 -> 0 (male)
            fsex=lambda df: 1 - df["msex"],
        )
        .drop(columns=["age_death", "msex"])
        .set_index("barcode")
    )

    # Save metadata for NB
    metadata.to_csv(out_path, sep="\t", index=True, compression="gzip")
    print(f"Saved metadata to: {out_path}")


if __name__ == "__main__":
    main()
