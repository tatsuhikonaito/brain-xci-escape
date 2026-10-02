# Bulk microglia (MiGA)

Sex differential expression and XCIR analysis across nine brain regions, using all ancestries. DGE uses dream with a participant random intercept; XCIR uses female samples.

- `01_DGE_sex.R`: Sex differential expression with SVA.
- `02_prepare_XCIR_input.sh`: Variant calling and population-specific 1KG filtering.
- `03_run_XCIR.R`: Sample-level XCIR beta-binomial models.
- `04_meta_XCIR.py`: Regional and participant-deduplicated XCIR meta-analysis.

DGE takes an all-gene count matrix, MiGA metadata and a sample key (`sample_id`, `participant_id`, `region`). XCIR meta-analysis takes metadata with `sample_id`, `participant_id`, `cohort`, `sex`, `region`, such as `MiGA_1.0_2.0_3.0_ALL.metadata_for_XCIR.txt`, and reads the sample-level results directly. Cohort labels are `MiGA_1.0`, `MiGA_2.0` and `MiGA_3.0`.

Requires R with dplyr, readr, edgeR, limma, variancePartition, sva, optparse, data.table and XCIR; samtools, bcftools and tabix; Python with numpy, pandas and scipy.

Run from the repository root; replace the example input paths:

```bash
Rscript bulk_microglia/01_DGE_sex.R counts.tsv metadata.tsv sample_key.tsv results/DGE ALL TRUE
bash bulk_microglia/02_prepare_XCIR_input.sh sample.bam reference.fa EUR.sites.vcf.gz processed_data/MiGA_1.0/sample
Rscript bulk_microglia/03_run_XCIR.R --vcf processed_data/MiGA_1.0/sample.chrX.XCIR_input.vcf.gz --outdir results/XCIR/MiGA_1.0/sample
python bulk_microglia/04_meta_XCIR.py --metadata metadata_for_XCIR.tsv --results-dir results/XCIR --vcf-dir processed_data --outdir results/XCIR/meta
```

Repeat the VCF and XCIR steps per sample, using the corresponding ancestry-specific sites VCF. The meta-analysis retains one sample per participant based on the largest number of XCIR-input variants, then combines `tau` by inverse-variance weighting. Escape requires BH FDR < 0.05, mean `tau` > 0.10 and at least three usable estimates. A metadata subset restricted to >25:75-skewed samples can be supplied for the sensitivity analysis.
