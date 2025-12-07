#!/usr/bin/env python3
"""
Prepare per-cell metadata for NB regression (MIT_ROSMAP only).

- Computes nUMI from raw counts
- Performs log-normalization, PCA, and Harmony (batch = individual_ID)
- Extracts oriPC_1-10 and hmPC_1-10
- Merges with clinical metadata
- Writes out a metadata table for downstream NB analyses

Environment variables
---------------------
XCI_BASE_DIR      : base directory for XCI project (contains <cohort>/...)
BIGBRAIN_BASE_DIR : base directory for snRNA-seq data (contains <cohort>/...)
(optional) XCI_REF_DIR : reference directory; if not set, defaults to
                         ${XCI_BASE_DIR}/ref

Output
------
${XCI_BASE_DIR}/MIT_ROSMAP/data/MIT_ROSMAP.metadata.txt.gz
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import harmonypy as hm
import scipy.sparse as sp


def main():
    # ------------------------------------------------
    # Settings (MIT_ROSMAP only)
    # ------------------------------------------------
    cohort = "MIT_ROSMAP"

    xci_base_dir = os.environ.get("XCI_BASE_DIR")
    bigbrain_base_dir = os.environ.get("BIGBRAIN_BASE_DIR")
    ref_dir = os.environ.get("XCI_REF_DIR")

    if xci_base_dir is None:
        raise RuntimeError("Please set XCI_BASE_DIR environment variable.")
    if bigbrain_base_dir is None:
        raise RuntimeError("Please set BIGBRAIN_BASE_DIR environment variable.")
    if ref_dir is None:
        ref_dir = os.path.join(xci_base_dir, "ref")

    base_dir = os.path.join(xci_base_dir, cohort)
    data_dir = os.path.join(bigbrain_base_dir, cohort)

    os.makedirs(os.path.join(base_dir, "data"), exist_ok=True)
    os.chdir(base_dir)

    # ------------------------------------------------
    # 1) Gene annotation (chrX) – kept for completeness
    # ------------------------------------------------
    annot_path = os.path.join(ref_dir, "genes.chrX.gtf.gz")
    if os.path.exists(annot_path):
        ref_annot = pd.read_csv(
            annot_path,
            sep="\t",
            header=None,
            names=[
                "seqname",
                "source",
                "feature",
                "start",
                "end",
                "score",
                "strand",
                "frame",
                "attribute",
            ],
        )
        ref_annot = ref_annot.loc[ref_annot["feature"] == "gene"].copy()

        # gene_name
        ref_annot["gene_name"] = [
            attr.split(";")[3]
            .replace("gene_name ", "")
            .replace('"', "")
            .lstrip(" ")
            for attr in ref_annot["attribute"]
        ]

        # index = gene_id (without version)
        ref_annot.index = [
            attr.split(";")[0]
            .replace("gene_id ", "")
            .replace('"', "")
            .split(".")[0]
            for attr in ref_annot["attribute"]
        ]
        ref_annot = ref_annot.loc[~ref_annot.index.duplicated()]
    else:
        # Not strictly needed for the rest of the script; safe to skip if absent
        ref_annot = None

    # ------------------------------------------------
    # 2) Read raw AnnData (MITROSMAP.h5ad)
    # ------------------------------------------------
    adata_path = os.path.join(
        data_dir, "analysis/snRNAseq/scanpy", "MITROSMAP.h5ad"
    )
    if not os.path.exists(adata_path):
        raise FileNotFoundError(f"h5ad file not found: {adata_path}")

    adata = sc.read_h5ad(adata_path)

    # ------------------------------------------------
    # 3) Total UMI from raw counts (for NB covariate)
    # ------------------------------------------------
    X = adata.X
    if sp.issparse(X):
        nUMI = np.asarray(X.sum(axis=1)).ravel()
    else:
        nUMI = X.sum(axis=1)

    # log10(nUMI + 1)
    adata.obs["nUMI"] = np.log10(nUMI + 1.0)

    # ------------------------------------------------
    # 4) Log-normalization for PCA / Harmony
    # ------------------------------------------------
    # library-size normalization
    sc.pp.normalize_total(adata, target_sum=1e4)
    # log-transform
    sc.pp.log1p(adata)

    # (Optional) HVG + scaling
    # sc.pp.highly_variable_genes(adata, n_top_genes=3000)
    # adata = adata[:, adata.var["highly_variable"]]
    # sc.pp.scale(adata, max_value=10)

    # ------------------------------------------------
    # 5) PCA & Harmony
    # ------------------------------------------------
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

    # ------------------------------------------------
    # 6) Build per-cell metadata (projid / individualID / nUMI / Annotation / PCs)
    # ------------------------------------------------
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

    # ------------------------------------------------
    # 7) Merge with clinical metadata
    # ------------------------------------------------
    clinical_path = os.path.join(base_dir, "data", "ROSMAP_clinical.csv")
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
    )

    # ------------------------------------------------
    # 8) Save metadata for NB
    # ------------------------------------------------
    out_path = os.path.join(base_dir, "data", f"{cohort}.metadata.txt.gz")
    metadata.to_csv(out_path, sep="\t", index=True, compression="gzip")
    print(f"Saved metadata to: {out_path}")


if __name__ == "__main__":
    main()
