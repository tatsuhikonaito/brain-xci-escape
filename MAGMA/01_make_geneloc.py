#!/usr/bin/env python3
import gzip
from pathlib import Path
import pandas as pd

# Input: GENCODE GTF (chrX)
GTF_PATH = Path("ref/gencode.v38.primary_assembly.annotation.gtf.gz")
REF_DIR = Path("/path/to/magma_ref")

# Output: MAGMA gene location file (chrX)
OUT_GENE_LOC = REF_DIR / "gencode.v38.chrX.MAGMA.gene.loc"


def parse_gene_name(attrs: str) -> str | None:
    """Extract gene_name from the 9th GTF column."""
    for item in attrs.split(";"):
        item = item.strip()
        if item.startswith("gene_name"):
            parts = item.split('"')
            if len(parts) >= 2:
                return parts[1]
            return None
    return None


def build_gene_df(gtf_path: Path) -> pd.DataFrame:
    """Read a GTF and return unique chrX genes as MAGMA-ready coordinates."""
    rows = []
    with gzip.open(gtf_path, "rt") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue

            chrom, _, feature, start, end, _, strand, _, attrs = fields
            if feature != "gene":
                continue
            if chrom not in {"chrX", "X"}:
                continue

            gene_name = parse_gene_name(attrs)
            if gene_name is None:
                continue

            rows.append((23, int(start), int(end), strand, gene_name))

    df = pd.DataFrame(rows, columns=["CHR", "START", "END", "STRAND", "SYMBOL"])

    # Keep one entry per gene symbol
    df = df.drop_duplicates(subset=["SYMBOL"], keep="first").sort_values(
        ["CHR", "START"]
    )
    df = df.reset_index(drop=True)
    return df


def write_gene_loc(df: pd.DataFrame, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)

    df_out = df.copy()
    df_out.insert(0, "MAGMA_ID", range(1, len(df_out) + 1))
    df_out = df_out[["MAGMA_ID", "CHR", "START", "END", "STRAND", "SYMBOL"]]

    # MAGMA gene.loc has no header
    df_out.to_csv(out_path, sep="\t", index=False, header=False)
    print(f"[DONE] wrote {len(df_out)} genes -> {out_path}")


def main() -> None:
    df = build_gene_df(GTF_PATH)
    write_gene_loc(df, OUT_GENE_LOC)


if __name__ == "__main__":
    main()
