# Single-cell differential expression

Negative-binomial mixed models for sex effects and sex-by-cell-state interactions.

- `01_prepare_NB_metadata.py`: Prepare UMI covariates, PCA, Harmony components and clinical metadata.
- `02_run_NB_sex_assoc.R`: Fit mixed models with a participant random intercept and perform likelihood-ratio tests.

Set input and output paths in both scripts. Use raw counts in `adata.X`, a cell-by-gene raw-count TSV with matching `barcode` values, and the participant list used for pseudobulk analysis.
Run from the repository root:

```bash
python scDGE/01_prepare_NB_metadata.py
Rscript scDGE/02_run_NB_sex_assoc.R
```

Requires Python with numpy, pandas, scipy, scanpy and harmonypy; R with lme4, Matrix and tidyverse.
