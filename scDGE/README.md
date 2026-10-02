# Single-cell differential expression

Negative-binomial mixed models for sex effects and sex-by-cell-state interactions, using MIT_ROSMAP as the example dataset. The models include a participant random intercept and use the first ten Harmony components to represent cell state.

- `01_prepare_NB_metadata.py`: Compute UMI covariates, PCA and Harmony components, then join clinical metadata.
- `02_run_NB_sex_assoc.R`: Fit sex-effect and interaction models and perform likelihood-ratio tests.

Requires Python with numpy, pandas, scipy, scanpy and harmonypy; R with lme4, Matrix and tidyverse. Metadata preparation expects raw counts in `adata.X`. Regression uses a cell-by-gene raw-count TSV with a `barcode` column and a participant list matching the pseudobulk analysis.

Run from the repository root:

```bash
export XCI_BASE_DIR=/path/to/XCI
export BIGBRAIN_BASE_DIR=/path/to/bigbrain/data
python scDGE/01_prepare_NB_metadata.py
# Set base_dir and the input paths in the R script.
Rscript scDGE/02_run_NB_sex_assoc.R
```

Metadata are saved under `MIT_ROSMAP/data/`; coefficients and likelihood-ratio test results are saved under `MIT_ROSMAP/NB/results/`. Barcodes must agree between the metadata and expression tables.
