# Bulk microglial RNA-seq (MiGA)

Representative scripts for sex differential expression and XCIR analysis. DGE and
XCIR are separate branches. Local file paths and cluster-specific job submission
are replaced by explicit inputs; the statistical calculations retain the study
settings.

| Script | Purpose |
| --- | --- |
| `01_DGE_sex.R` | Sex differential expression using dream across nine brain regions |
| `02_prepare_XCIR_input.sh` | ChrX variant calling and two-tier 1KG-based filtering |
| `03_run_XCIR.R` | Single-sample XCIR beta-binomial models |
| `04_collect_XCIR_results.py` | Collection of explicitly selected sample results by region |
| `05_integrate_XCIR.py` | Per-region and donor-deduplicated XCIR meta-analysis |

Dependencies: R with dplyr, readr, tibble, edgeR, limma, variancePartition, sva,
optparse, data.table and XCIR; samtools, bcftools and tabix; Python with numpy,
pandas and scipy. Use the original software versions and reference files for
numerical reproduction. Run information is saved with the outputs. XCIR annotation
may require network access.

**Run all commands from the repository root, where `bulk_microglia/` is visible.**
All paths, including paths inside TSV manifests, refer to that working directory
unless absolute. The `bulk_microglia/data/` filenames below are examples, not bundled data.

## 1. Sex differential expression

Inputs are a gene-by-sample count TSV (first column: gene IDs), a metadata TSV, and
a sample-key TSV (`sample_id`, `participant_id`, `region`; one row per sample).
Normalized regions are `MFG`, `STG`, `SVZ`, `THA`, `CC`, `CER`, `HIP`, `OCC`, `SN`.
Supply the original complete expression matrix, not an X-only subset.

Metadata columns: `sample_id`, `Sex` (`Male`/`Female`), `Age`, `Ancestry`,
`Diagnosis`, `APOE`, `Cause_of_Death`, `PMD_min`, `PCT_MRNA_BASES`,
`PCT_RIBOSOMAL_BASES`, `MEDIAN_5PRIME_TO_3PRIME_BIAS`, `Cohort`. The key supplies
`participant_id` and `region`; `Cohort` comes from metadata, not a key filename.

```bash
Rscript bulk_microglia/01_DGE_sex.R bulk_microglia/data/MiGA.counts.tsv bulk_microglia/data/MiGA.metadata.tsv bulk_microglia/data/MiGA.sample_key.tsv bulk_microglia/results/MiGA/DGE ALL TRUE
```

Defaults are all ancestries and SVA. The model retains multiple samples per donor,
with a participant random intercept and region as a fixed effect. Positive
`SexFemale` logFC indicates higher female expression. Count rounding,
`filterByExpr`, TMM renormalization, scaling of five numeric covariates, voom-based
SVA (`num.sv` with `leek`, `sva`, null model `~1`), and
`voomWithDreamWeights` / `dream` / `eBayes` retain their original order.

Diagnostic/cause-of-death recoding and complete-case selection are preserved.
`APOE4` is 1 for 24/34/44, 0 for 23/33, and missing otherwise. The `Ancestry` term
is not automatically dropped in ancestry-restricted runs; verify the design is
estimable before using such a setting. Outputs include the original-style DGE
filename, retained design table and session information. No gene-ID remapping or
additional X/non-PAR filtering is applied.

## 2. Single-sample XCIR

Supply an indexed BAM, its matching indexed reference FASTA, and the original
population-specific 1KG sites VCF. Use consistent contig naming in all three.
Reference-panel preparation and ancestry assignment are upstream inputs.

```bash
bash bulk_microglia/02_prepare_XCIR_input.sh bulk_microglia/data/bams/sample01.bam bulk_microglia/data/reference.fa bulk_microglia/data/1KG/EUR.sites.vcf.gz bulk_microglia/processed_data/MiGA/sample01
Rscript bulk_microglia/03_run_XCIR.R --vcf bulk_microglia/processed_data/MiGA/sample01.chrX.XCIR_input.vcf.gz --outdir bulk_microglia/results/MiGA/XCIR/sample01 --read_count_cutoff 6 --het_cutoff 1 --release hg38 --model BB
```

For an `X` BAM, use the generated `.X.XCIR_input.vcf.gz` path. For AFR/AMR/EAS,
provide the corresponding sites VCF instead of EUR; do not infer an unknown
ancestry. Prepare the shared sites index before launching concurrent jobs.

Retained settings: MAPQ >=20, base quality >=20; Tier 1 DP >=6 and both allele
counts >=1; Tier 2 DP >=20 and both allele counts >=3; VCF DP <=100000; biallelic
SNPs. Variants are called first, partitioned by the sites panel, filtered and
concatenated. No extra genotype or PAR filter is added. `MAX_DP=100000` is a VCF
filter, **not** an mpileup `-d` setting; the original mpileup defaults are retained.

XCIR uses `readVCF4`, `annotateX`, `getGenicDP(highest_expr=TRUE)` and
`betaBinomXI(model="BB")`. The operative annotation cutoffs are 6/1. Intermediate
tables, the separate `sample_clean` output and run information are saved. Each
invocation recomputes results rather than silently reusing an existing fit.
`getXCIstate` is not used in the gene-level meta-analysis.

## 3. Collect analysis-ready results

**The analysis sample list is an explicit input.** Collection and integration do
not infer sex/skewness eligibility, gene-ID harmonization, or a `sample_clean`
filter. Reproduce these upstream choices before assembling the input tables.

With original regional `betaBinomXI` tables available, skip collection and use
them directly in step 4. Otherwise prepare a TSV with `sample_id`, `region`,
`result_file`, listing exactly the analysis-eligible sample results. The `sample`
column inside each result must match `sample_id`. File 04 is a mechanical adapter:
it concatenates the specified results without further selection or ID conversion.

```bash
python bulk_microglia/04_collect_XCIR_results.py --manifest bulk_microglia/data/MiGA.xcir_samples.tsv --outdir bulk_microglia/results/MiGA/XCIR/regional
```

## 4. Gene-level integration

`--region-files` is a TSV with `region`, `path`, one row per regional result table.
Use the order `MFG`, `STG`, `SVZ`, `THA`, `CC`, `CER`, `HIP`, `OCC`, `SN` for the
extended-region analysis. Each table needs `GENE`, `sample`, `tau`, `var_tau`.
The collector writes `region_files.tsv` in first-appearance order; alternatively
prepare it to point to the original regional results.

`--sample-key` needs `sample_id`, `participant_id`, `vcf`, with one row per sample
and its exact XCIR-input VCF path. Extra columns are allowed. It may also serve as
the DGE key when it includes `region`. Every analyzed sample must have a mapping;
extra key rows do not add samples to the analysis.

```bash
python bulk_microglia/05_integrate_XCIR.py --region-files bulk_microglia/results/MiGA/XCIR/regional/region_files.tsv --sample-key bulk_microglia/data/MiGA.xcir_sample_key.tsv --outdir bulk_microglia/results/MiGA/XCIR/meta --prefix MiGA.ALL
```

Per-region summaries retain all supplied regional result rows. The combined
summary selects **one sample per participant**, maximizing the number of records
in the whole XCIR-input VCF among samples represented in the results. This is not
a gene-specific selection. Ties retain the original input-order-dependent behavior;
preserve both region order and within-region row order. Selected samples are saved.

For each gene, inverse-variance weighting uses `tau` and positive, non-missing
`var_tau`. Two-sided normal P values are BH-adjusted **only among genes with
`n_meta >= 3`**. `rho_mean` is the unweighted mean of non-missing **raw `tau`**,
including estimates without a valid variance; no clipping or separate `n_rho`
threshold is applied. The XCIR call requires `padj_meta < 0.05`, `rho_mean > 0.10`
and `n_meta >= 3`.

Outputs: `{prefix}.{region}.XCIR_gene_meta.tsv.gz`,
`{prefix}.ALLREGIONS.XCIR_gene_meta.tsv.gz`, `{prefix}.selected_samples.tsv` and
`{prefix}.run_info.json`. Column names follow the study integration code. Missing
files/mappings and malformed numeric strings cause errors; genuine missing
estimates retain the original handling.

These scripts produce DGE and XCIR statistics, not the final annotated manuscript
table or a new joint DGE-plus-XCIR classification. Exact study reproduction also
requires matching sample membership, references, identifiers, row order and software
versions; compare the original regional inputs and selected-sample list first.
