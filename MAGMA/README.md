# MAGMA gene-set analysis

Gene and gene-set tests for X-chromosomal AD association signals.

- `01_make_geneloc.py`: Create gene coordinates and a symbol-to-ID mapping from GENCODE v38.
- `02_make_geneset.py`: Define escape/non-escape gene sets, including separate bulk microglia sets.
- `03_run_MAGMA.sh`: Prepare summary statistics and run MAGMA.

Requires Python with pandas and MAGMA v1.10. Use a prepared female EUR non-PAR LD reference in GRCh38.
Set input paths and summary-statistic columns in the scripts, then run from the repository root:

```bash
python MAGMA/01_make_geneloc.py
python MAGMA/02_make_geneset.py
bash MAGMA/03_run_MAGMA.sh GCST90444373
bash MAGMA/03_run_MAGMA.sh GCST90449045
```

Outputs are written to `MAGMA/data/` and `MAGMA/results/` under `BASE_DIR`.
