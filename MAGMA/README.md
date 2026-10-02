# MAGMA gene-set analysis

Gene and gene-set tests for X-chromosomal AD association signals.

- `01_make_geneloc.py`: Create gene coordinates and a symbol-to-ID mapping from GENCODE v38.
- `02_make_geneset.py`: Define cell-type gene sets from Xi-ratio estimates and previous XCI annotations.
- `03_run_MAGMA.sh`: Prepare summary statistics and run MAGMA.

The gene sets use Xi ratio > 0.1 **or** a previous `nonPAR escape` annotation, excluding XIST; DGE significance is not required. The summary table and all reference files must be prepared beforehand. The LD reference is a female European, non-PAR X-chromosome PLINK panel in GRCh38.

Requires Python with pandas and MAGMA. Set paths in the three scripts, then run from the repository root:

```bash
python MAGMA/01_make_geneloc.py
python MAGMA/02_make_geneset.py
bash MAGMA/03_run_MAGMA.sh GCST90444373
bash MAGMA/03_run_MAGMA.sh GCST90449045
```

Summary-statistic column positions and `HAS_HEADER` are set near the top of `03_run_MAGMA.sh`; adjust these to the input layout. Gene-set IDs must use the same gene-location file as the MAGMA annotation step. Results are written to `MAGMA/results/` under `BASE_DIR`.
