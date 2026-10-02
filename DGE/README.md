# Pseudobulk differential expression

`pseudobulk_DGE.R` fits DESeq2 models for sex differences or AD-related expression in MIT_ROSMAP, ROSMAP and SEA_AD. SVA is enabled by default. The `pheno` analysis runs separately in females, males and both sexes.

Requires R with dplyr, readr, DESeq2 and sva. Inputs are pseudobulk RData files containing `genes_counts` (first column: gene IDs; remaining columns: sample counts) and cohort-specific clinical metadata. For `celltype=all`, the combined pseudobulk file must already exist.

Run from the repository root, using absolute paths for the data directories:

```bash
export XCI_BASE_DIR=/path/to/XCI
export BIGBRAIN_BASE_DIR=/path/to/bigbrain/data
Rscript DGE/pseudobulk_DGE.R MIT_ROSMAP MG sex AD_NCI TRUE
Rscript DGE/pseudobulk_DGE.R MIT_ROSMAP MG pheno AD_NCI TRUE
```

Arguments: `cohort celltype analysis pheno_definition [adjust_batch]`. Phenotype definitions are `AD_else`, `AD_NCI` or `ADMCI_NCI` (the last is ROSMAP-based only); this argument is ignored for `sex`.

ROSMAP also uses `data/sample_list.ROSMAP.notin_MIT_ROSMAP.tsv` to run an analysis excluding overlapping participants. Results are written to `<XCI_BASE_DIR>/<cohort>/DEG/results/`. Input paths and cohort-specific coding are defined in the script.
