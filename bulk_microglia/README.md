# Bulk microglia

Scripts for the bulk microglial RNA-seq analyses in MiGA.

- `01_DGE_sex.R`: sex differential expression using `dream`.
- `02_prepare_XCIR_input.sh`: prepare chrX variants for XCIR.
- `03_run_XCIR.R`: run XCIR for each sample.
- `04_meta_XCIR.py`: gene-level XCIR meta-analysis across samples and regions.

Input paths are supplied to each script. Analysis settings are defined in the code.
