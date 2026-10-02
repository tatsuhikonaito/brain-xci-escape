#!/usr/bin/env Rscript

# MiGA bulk microglia: sex differential expression across nine brain regions.
# Keep all eligible regional samples; model repeated donors with dream.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(edgeR)
  library(limma)
  library(variancePartition)
  library(sva)
})

# Run from the repository root; paths are relative to the working directory.
# The prepared key has sample_id, participant_id and normalized region columns.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4 || length(args) > 6) {
  stop("Usage: Rscript bulk_microglia/01_DGE_sex.R <counts.tsv> <metadata.tsv> <sample_key.tsv> <outdir> [ancestry=ALL] [adjust_batch=TRUE]")
}
genes_counts_path <- args[1]
meta_path <- args[2]
key_path <- args[3]
dir_out <- args[4]
ancestry <- if (length(args) >= 5) args[5] else "ALL"
adjust_batch <- if (length(args) >= 6) as.logical(args[6]) else TRUE
if (is.na(adjust_batch)) stop("adjust_batch must be TRUE or FALSE")
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

keys <- read_tsv(key_path, show_col_types = FALSE) %>%
  transmute(sample_id = as.character(sample_id),
            participant_id = as.character(participant_id),
            region = as.character(region)) %>%
  distinct()
if (anyDuplicated(keys$sample_id)) stop("Conflicting sample_id rows in sample_key")

metadata <- read_tsv(meta_path, show_col_types = FALSE)

metadata <- metadata %>%
  mutate(
    Dx_group = case_when(
      grepl("Control|Healthy|Non[_ -]?Demented[_ -]?Control", Diagnosis, ignore.case = TRUE) ~ "Control",
      grepl("AD|Alzheimer|Dementia|LBD|Lewy|PD|Parkinson|PSP|FTD|MCI|Neurodegeneration|Multi_Atrophy|Vascular dementia",
            Diagnosis, ignore.case = TRUE) ~ "Neurodegenerative",
      grepl("Depression|Bipolar|Schizo|Psychotic", Diagnosis, ignore.case = TRUE) ~ "Psychiatric",
      TRUE ~ "Other"
    ),
    Dx_group = factor(Dx_group, levels = c("Control", "Psychiatric", "Neurodegenerative", "Other")),
    APOE4 = ifelse(APOE %in% c("24","34","44"), 1, ifelse(APOE %in% c("23","33"), 0, NA)),
    COD_group = case_when(
      grepl("Cancer", Cause_of_Death, ignore.case = TRUE) ~ "Cancer",
      grepl("Cardio|Cardiorespiratory|heart|resp", Cause_of_Death, ignore.case = TRUE) ~ "Cardiorespiratory",
      grepl("Infect", Cause_of_Death, ignore.case = TRUE) ~ "Infection",
      grepl("Neurodegeneration", Cause_of_Death, ignore.case = TRUE) ~ "Neurodegeneration",
      grepl("Suicide", Cause_of_Death, ignore.case = TRUE) ~ "Suicide",
      TRUE ~ "Other"
    ),
    COD_group = factor(COD_group),
    Sex = relevel(factor(Sex), ref = "Male"),
    Ancestry = factor(Ancestry)
  )

genes_counts <- read.table(
  genes_counts_path,
  header = TRUE, sep = "\t", row.names = 1, check.names = FALSE
)

# The key defines region and participant_id; Cohort comes from metadata.
data0 <- keys %>%
  inner_join(metadata, by = "sample_id") %>%
  filter(Ancestry == ancestry | ancestry == "ALL")

samples_use <- intersect(colnames(genes_counts), data0$sample_id)
data0 <- data0 %>% filter(sample_id %in% samples_use) %>% arrange(match(sample_id, samples_use))
counts0 <- genes_counts[, samples_use, drop = FALSE]

covariates <- c(
  "Age",
  "Ancestry",
  "Dx_group",
  "APOE4",
  "COD_group",
  "PMD_min",
  "PCT_MRNA_BASES",
  "PCT_RIBOSOMAL_BASES",
  "MEDIAN_5PRIME_TO_3PRIME_BIAS",
  "Cohort",
  "region"
)

vars_needed <- c("Sex", "participant_id", covariates)
keep <- apply(!is.na(data0[, vars_needed, drop = FALSE]), 1, all)
data <- data0[keep, , drop = FALSE]
counts <- counts0[, keep, drop = FALSE]

data$Cohort <- factor(data$Cohort)
data$region <- factor(data$region, levels = c("MFG", "STG", "SVZ", "THA", "CC", "CER", "HIP", "OCC", "SN"))
data$participant_id <- factor(data$participant_id)

# Count rounding, normalization and expression filtering match the study code.
dge <- DGEList(counts = round(counts))
dge <- calcNormFactors(dge)

fixed_terms <- paste(c("Sex", covariates), collapse = " + ")
design_fixed <- model.matrix(as.formula(paste("~", fixed_terms)), data = data)

keep_gene <- filterByExpr(dge, design_fixed)
dge <- dge[keep_gene, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)

# Align metadata rows to expression columns (must be consistent for dream)
data <- data[match(colnames(dge$counts), data$sample_id), , drop = FALSE]
rownames(data) <- data$sample_id

# Rescale numeric covariates to avoid "different scales" model failures
num_vars <- c("Age", "PMD_min", "PCT_MRNA_BASES", "PCT_RIBOSOMAL_BASES", "MEDIAN_5PRIME_TO_3PRIME_BIAS")
num_vars <- intersect(num_vars, colnames(data))
data[, num_vars] <- lapply(data[, num_vars, drop = FALSE], function(x) as.numeric(scale(as.numeric(x))))

# Estimate SVs from voom E using a fixed-effect model, then include them as covariates
sv_terms <- character(0)
v_fixed <- voom(dge, design_fixed, plot = FALSE)

if (adjust_batch) {
  mod <- design_fixed
  mod0 <- model.matrix(~ 1, data = data)
  nsv <- num.sv(v_fixed$E, mod, method = "leek")
  if (!is.na(nsv) && nsv > 0) {
    sv <- sva(v_fixed$E, mod, mod0, n.sv = nsv)$sv
    colnames(sv) <- paste0("SV", seq_len(ncol(sv)))
    data <- cbind(data, sv)
    sv_terms <- colnames(sv)
  }
}

# Fit dream model with random intercept for repeated measures (participant_id)
rhs <- c("Sex", covariates)
if (length(sv_terms) > 0) rhs <- c(sv_terms, rhs)
form <- as.formula(paste("~", paste(rhs, collapse = " + "), "+ (1|participant_id)"))

vobj <- voomWithDreamWeights(dge, form, data)
fit <- eBayes(dream(vobj, form, data))

base_out <- paste0("MiGA.", ancestry, ".ALL.extended_regions.DEG_sex.dream", ifelse(adjust_batch, ".sva", ""))
outfile_main <- file.path(dir_out, paste0(base_out, ".txt.gz"))

res_main <- topTable(fit, coef = "SexFemale", number = Inf, sort.by = "P") %>%
  tibble::rownames_to_column("gene")

write.table(res_main, gzfile(outfile_main), quote = FALSE, sep = "\t", row.names = FALSE)


# Save the retained design and software versions without changing the fit.
write.table(data, gzfile(file.path(dir_out, paste0(base_out, ".design.tsv.gz"))),
            quote = FALSE, sep = "\t", row.names = FALSE)
capture.output(sessionInfo(), file = file.path(dir_out, paste0(base_out, ".sessionInfo.txt")))
