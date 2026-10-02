#!/usr/bin/env Rscript

# Single-sample XCIR analysis.
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(XCIR)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--vcf", type = "character", help = "Filtered chrX VCF"),
  make_option("--outdir", type = "character", help = "Sample output directory"),
  make_option("--read_count_cutoff", type = "integer", default = 6),
  make_option("--het_cutoff", type = "integer", default = 1),
  make_option("--release", type = "character", default = "hg38"),
  make_option("--model", type = "character", default = "BB")
)))
if (is.null(opt$vcf) || is.null(opt$outdir)) stop("--vcf and --outdir are required")
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

write_tab <- function(dt, name) {
  con <- gzfile(file.path(opt$outdir, name), "wt")
  on.exit(close(con))
  write.table(dt, con, sep = "\t", quote = FALSE, row.names = FALSE)
}

snp_dt <- readVCF4(opt$vcf)
if ("CHROM" %in% names(snp_dt)) snp_dt[, CHROM := gsub("^chr", "", CHROM)]

anno_dt <- annotateX(snp_dt, read_count_cutoff = opt$read_count_cutoff,
                     het_cutoff = opt$het_cutoff, release = opt$release)

genic_dt <- getGenicDP(anno_dt, highest_expr = TRUE)

bb_dt <- betaBinomXI(genic_dt, model = opt$model)
write_tab(bb_dt, "betaBinomXI_results.tsv.gz")
