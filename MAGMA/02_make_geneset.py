#!/usr/bin/env python3
import os
from pathlib import Path
import pandas as pd

# =========================================
# Paths / settings (match gene.loc script style)
# =========================================
BASE_DIR = Path("/path/to/XCI_project")
os.chdir(BASE_DIR)

SUMMARY_PATH = Path("/path/to/summary")

REF_DIR = Path("/path/to/magma_ref")
OUT_GENE_LOC = REF_DIR / "gencode.v38.chrX.MAGMA.gene.loc"

# Output: MAGMA gene-set file (.sets)
OUT_GENESETS = BASE_DIR / "MAGMA/data/magma_genesets.txt"

# Columns in summary_df
GENE_COL = "Gene"
CELLTYPE_COL = "Cell type"
STATUS_COL = "Annotation of XCI status"
FDR_META_COL = "FDR (ROSMAP_MIT_ROSMAP)"
XI_META_COL = "Ratio of the expression from Xi (ROSMAP_MIT_ROSMAP)"

CELLTYPES = ["Ast", "Exc", "Inh", "Mic", "Oli", "OPC"]

FDR_THR = 0.05
XI_THR = 0.1


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
    df["deg_sig"] = df[FDR_META_COL].notna() & (df[FDR_META_COL] < FDR_THR)
    df["xi_measured"] = df[XI_META_COL].notna()
    df["xi_escape"] = df["xi_measured"] & (df[XI_META_COL] > XI_THR)
    return df


def genes_to_magma_ids(gene_symbols: set[str], symbol_to_id: dict[str, int]) -> list[int]:
    ids = sorted({symbol_to_id[g] for g in gene_symbols if g in symbol_to_id})
    return ids


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

    # ----- per-cell-type: xi_or_annotEscape_{ct} -----
    for ct in CELLTYPES:
        sub = df[df[CELLTYPE_COL] == ct].copy()
        if sub.empty:
            continue

        genes_annot_escape = set(sub.loc[sub[STATUS_COL] == "nonPAR escape", GENE_COL].dropna().astype(str))
        genes_xi_escape = set(sub.loc[sub["xi_escape"], GENE_COL].dropna().astype(str))

        genes_xi_or_annot = genes_xi_escape.union(genes_annot_escape)
        ids = genes_to_magma_ids(genes_xi_or_annot, symbol_to_id)

        if ids:
            genesets[f"xi_or_annotEscape_{ct}"] = ids

    # ----- anyCT: xi_or_annotEscape_anyCT -----
    genes_annot_escape_any = set(df.loc[df[STATUS_COL] == "nonPAR escape", GENE_COL].dropna().astype(str))
    genes_xi_escape_any = set(df.loc[df["xi_escape"], GENE_COL].dropna().astype(str))

    genes_xi_or_annot_any = genes_xi_escape_any.union(genes_annot_escape_any)
    ids_any = genes_to_magma_ids(genes_xi_or_annot_any, symbol_to_id)
    if ids_any:
        genesets["xi_or_annotEscape_anyCT"] = ids_any

    # Write MAGMA .sets (no header)
    OUT_GENESETS.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_GENESETS, "w") as out:
        for set_name in sorted(genesets.keys()):
            out.write(set_name + "\t" + "\t".join(map(str, genesets[set_name])) + "\n")

    print(f"[DONE] wrote {len(genesets)} gene sets -> {OUT_GENESETS}")


if __name__ == "__main__":
    main()
