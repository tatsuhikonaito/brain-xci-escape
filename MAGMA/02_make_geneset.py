#!/usr/bin/env python3
import os
from pathlib import Path
import pandas as pd

# =========================================
# Paths / settings
# =========================================
BASE_DIR = Path("/path/to/XCI_project")
os.chdir(BASE_DIR)

SUMMARY_PATH = Path("/path/to/summary")
MICROGLIA_SUMMARY_PATH = Path("/path/to/summary_df.microglia.txt")

REF_DIR = Path("/path/to/magma_ref")
OUT_GENE_LOC = REF_DIR / "gencode.v38.chrX.MAGMA.gene.loc"

# Output: MAGMA gene-set file (.sets)
OUT_GENESETS = BASE_DIR / "MAGMA/data/magma_genesets.txt"
OUT_MICROGLIA_GENESETS = BASE_DIR / "MAGMA/data/magma_genesets.microglia_xcir_only.sets"

# Columns in summary_df
GENE_COL = "Gene"
CELLTYPE_COL = "Cell type"
STATUS_COL = "Annotation of XCI status"
XI_META_COL = "Ratio of the expression from Xi (ROSMAP_MIT_ROSMAP)"

CELLTYPES = ["Ast", "Exc", "Inh", "Mic", "Oli", "OPC"]

XI_THR = 0.1
XCIR_PADJ_THR = 0.05


# =========================================
# Helpers
# =========================================
def load_symbol_to_id(gene_loc_path: Path) -> dict[str, int]:
    gene_loc = pd.read_csv(gene_loc_path, sep="\t", header=None)
    gene_loc.columns = ["MAGMA_ID", "CHR", "START", "END", "STRAND", "SYMBOL"]

    symbol_to_id = (
        gene_loc.dropna(subset=["SYMBOL"])
        .drop_duplicates(subset=["SYMBOL"], keep="first")
        .set_index("SYMBOL")["MAGMA_ID"]
        .astype(int)
        .to_dict()
    )
    return symbol_to_id


def add_evidence_flags(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["xi_measured"] = df[XI_META_COL].notna()
    df["xi_escape"] = df["xi_measured"] & (df[XI_META_COL] > XI_THR)
    return df


def genes_to_magma_ids(gene_symbols: set[str], symbol_to_id: dict[str, int]) -> list[int]:
    ids = sorted({symbol_to_id[g] for g in gene_symbols if g in symbol_to_id})
    return ids


def write_genesets(genesets: dict[str, list[int]], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as out:
        for set_name in sorted(genesets):
            out.write(set_name + "\t" + "\t".join(map(str, genesets[set_name])) + "\n")
    print(f"[DONE] wrote {len(genesets)} gene sets -> {out_path}")


# =========================================
# Main
# =========================================
def main() -> None:
    # Load summary_df (drop XIST)
    df = pd.read_csv(SUMMARY_PATH, sep="\t")
    df = df[df[GENE_COL] != "XIST"].copy()

    # Load symbol -> MAGMA_ID mapping (from the same gene.loc you build elsewhere)
    symbol_to_id = load_symbol_to_id(OUT_GENE_LOC)

    # Evidence flags
    df = add_evidence_flags(df)

    genesets: dict[str, list[int]] = {}

    # ----- per-cell-type escape and non-escape sets -----
    for ct in CELLTYPES:
        sub = df[df[CELLTYPE_COL] == ct].copy()
        if sub.empty:
            continue

        genes_annot_escape = set(sub.loc[sub[STATUS_COL] == "nonPAR escape", GENE_COL].dropna().astype(str))
        genes_xi_escape = set(sub.loc[sub["xi_escape"], GENE_COL].dropna().astype(str))

        genes_xi_or_annot = genes_xi_escape.union(genes_annot_escape)
        ct_all_mappable = {g for g in sub[GENE_COL].unique() if g in symbol_to_id}
        genes_not_xi_or_annot = ct_all_mappable.difference(genes_xi_or_annot)

        for name, genes in [
            (f"xi_or_annotEscape_{ct}", genes_xi_or_annot),
            (f"not_xi_or_annotEscape_{ct}", genes_not_xi_or_annot),
        ]:
            ids = genes_to_magma_ids(genes, symbol_to_id)
            if ids:
                genesets[name] = ids

    write_genesets(genesets, OUT_GENESETS)

    # ----- bulk microglia, analyzed separately from the snRNA-seq Mic sets -----
    mic = pd.read_csv(MICROGLIA_SUMMARY_PATH, sep="\t")
    mic = mic[mic[GENE_COL] != "XIST"].copy()
    mic[GENE_COL] = mic[GENE_COL].astype("string")
    mic = mic[mic[GENE_COL].notna()].copy()
    for col in ["XCIR_padj_meta", "XCIR_rho_mean"]:
        mic[col] = pd.to_numeric(mic[col], errors="coerce")

    mic["xi_escape_mic"] = (
        mic["XCIR_padj_meta"].notna()
        & mic["XCIR_rho_mean"].notna()
        & (mic["XCIR_padj_meta"] < XCIR_PADJ_THR)
        & (mic["XCIR_rho_mean"] > XI_THR)
    )
    mic["annot_escape"] = mic["XCI_status"].astype("string").eq("nonPAR_escape")

    # Retain any-row support when a gene has multiple rows, as in the source analysis.
    gene_flags = (
        mic.groupby(GENE_COL, dropna=True)[["xi_escape_mic", "annot_escape"]]
        .max()
        .reset_index()
    )
    all_genes_mic = set(gene_flags[GENE_COL].astype(str))
    genes_xi_mic = set(gene_flags.loc[gene_flags["xi_escape_mic"], GENE_COL].astype(str))
    genes_annot_mic = set(gene_flags.loc[gene_flags["annot_escape"], GENE_COL].astype(str))
    genes_escape_mic = genes_xi_mic.union(genes_annot_mic)

    mic_genesets: dict[str, list[int]] = {}
    for name, genes in [
        ("xi_or_annotEscape_Mic", genes_escape_mic),
        ("not_xi_or_annotEscape_Mic", all_genes_mic.difference(genes_escape_mic)),
    ]:
        ids = genes_to_magma_ids(genes, symbol_to_id)
        if ids:
            mic_genesets[name] = ids

    write_genesets(mic_genesets, OUT_MICROGLIA_GENESETS)


if __name__ == "__main__":
    main()
