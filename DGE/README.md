# Pseudobulk differential expression

DESeq2 analyses of sex differences and AD-associated expression, with SVA by default.

Inputs: an RData file containing `genes_counts` (gene IDs in column 1, then sample counts) and a clinical CSV.
Set `clinical_format` in the script to `cogdx` or `diagnosis`, matching the metadata columns.
For an additional analysis excluding overlapping participants, set `nonoverlapping_samples_file` to a list of retained participant IDs.

Run from the repository root, using separate output prefixes for each dataset and cell type:

```bash
Rscript DGE/pseudobulk_DGE.R counts.RData clinical.csv results/Ast sex
Rscript DGE/pseudobulk_DGE.R counts.RData clinical.csv results/Ast pheno AD_NCI TRUE
```

Arguments: `counts_file clinical_file outprefix analysis [pheno_definition] [adjust_batch]`.
Requires R with dplyr, readr, DESeq2 and sva.
