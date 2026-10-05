#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(XCIR)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--vcf", type = "character"),
  make_option("--outdir", type = "character"),
  make_option("--read_count_cutoff", type = "integer", default = 6),
  make_option("--het_cutoff", type = "integer", default = 1),
  make_option("--release", type = "character", default = "hg38"),
  make_option("--model", type = "character", default = "BB")
)))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

write_tab <- function(dt, file) {
  con <- gzfile(file, "wt")
  on.exit(close(con))
  write.table(dt, con, sep = "\t", quote = FALSE, row.names = FALSE)
}

snp_dt <- readVCF4(opt$vcf)
if ("CHROM" %in% names(snp_dt)) snp_dt[, CHROM := gsub("^chr", "", CHROM)]

anno_dt <- annotateX(
  snp_dt,
  read_count_cutoff = opt$read_count_cutoff,
  het_cutoff = opt$het_cutoff,
  release = opt$release
)

genic_dt <- getGenicDP(anno_dt, highest_expr = TRUE)
bb_dt <- betaBinomXI(genic_dt, model = opt$model)
samp_dt <- sample_clean(bb_dt)

write_tab(bb_dt, file.path(opt$outdir, "betaBinomXI_results.tsv.gz"))
write_tab(samp_dt, file.path(opt$outdir, "sample_clean.tsv.gz"))
