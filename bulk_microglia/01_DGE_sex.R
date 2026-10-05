#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(edgeR)
  library(limma)
  library(variancePartition)
  library(sva)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript 01_DGE_sex.R <counts.tsv> <metadata.tsv> <sample_key.tsv> <outdir>")
}

counts_path <- args[1]
metadata_path <- args[2]
key_path <- args[3]
outdir <- args[4]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

keys <- read_tsv(key_path, show_col_types = FALSE) %>%
  transmute(
    sample_id = as.character(sample_id),
    participant_id = as.character(participant_id),
    region = as.character(region)
  ) %>%
  distinct()

metadata <- read_tsv(metadata_path, show_col_types = FALSE) %>%
  mutate(
    Dx_group = case_when(
      grepl("Control|Healthy|Non[_ -]?Demented[_ -]?Control", Diagnosis, ignore.case = TRUE) ~ "Control",
      grepl("AD|Alzheimer|Dementia|LBD|Lewy|PD|Parkinson|PSP|FTD|MCI|Neurodegeneration|Multi_Atrophy|Vascular dementia",
            Diagnosis, ignore.case = TRUE) ~ "Neurodegenerative",
      grepl("Depression|Bipolar|Schizo|Psychotic", Diagnosis, ignore.case = TRUE) ~ "Psychiatric",
      TRUE ~ "Other"
    ),
    Dx_group = factor(Dx_group, levels = c("Control", "Psychiatric", "Neurodegenerative", "Other")),
    APOE4 = ifelse(APOE %in% c("24", "34", "44"), 1, ifelse(APOE %in% c("23", "33"), 0, NA)),
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

counts_all <- read.table(
  counts_path,
  header = TRUE, sep = "\t", row.names = 1, check.names = FALSE
)

data0 <- keys %>% inner_join(metadata, by = "sample_id")
samples_use <- intersect(colnames(counts_all), data0$sample_id)
data0 <- data0 %>%
  filter(sample_id %in% samples_use) %>%
  arrange(match(sample_id, samples_use))
counts0 <- counts_all[, samples_use, drop = FALSE]

covariates <- c(
  "Age", "Ancestry", "Dx_group", "APOE4", "COD_group", "PMD_min",
  "PCT_MRNA_BASES", "PCT_RIBOSOMAL_BASES", "MEDIAN_5PRIME_TO_3PRIME_BIAS",
  "Cohort", "region"
)

vars_needed <- c("Sex", "participant_id", covariates)
keep <- apply(!is.na(data0[, vars_needed, drop = FALSE]), 1, all)
data <- data0[keep, , drop = FALSE]
counts <- counts0[, keep, drop = FALSE]

data$Cohort <- factor(data$Cohort)
data$region <- factor(data$region, levels = c("MFG", "STG", "SVZ", "THA", "CC", "CER", "HIP", "OCC", "SN"))
data$participant_id <- factor(data$participant_id)

dge <- DGEList(counts = round(counts))
dge <- calcNormFactors(dge)

fixed_terms <- paste(c("Sex", covariates), collapse = " + ")
design_fixed <- model.matrix(as.formula(paste("~", fixed_terms)), data = data)

keep_gene <- filterByExpr(dge, design_fixed)
dge <- dge[keep_gene, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)

data <- data[match(colnames(dge$counts), data$sample_id), , drop = FALSE]
rownames(data) <- data$sample_id

num_vars <- c("Age", "PMD_min", "PCT_MRNA_BASES", "PCT_RIBOSOMAL_BASES", "MEDIAN_5PRIME_TO_3PRIME_BIAS")
data[, num_vars] <- lapply(data[, num_vars, drop = FALSE], function(x) as.numeric(scale(as.numeric(x))))

v_fixed <- voom(dge, design_fixed, plot = FALSE)
mod <- design_fixed
mod0 <- model.matrix(~ 1, data = data)
nsv <- num.sv(v_fixed$E, mod, method = "leek")

sv_terms <- character(0)
if (!is.na(nsv) && nsv > 0) {
  sv <- sva(v_fixed$E, mod, mod0, n.sv = nsv)$sv
  colnames(sv) <- paste0("SV", seq_len(ncol(sv)))
  data <- cbind(data, sv)
  sv_terms <- colnames(sv)
}

rhs <- c(sv_terms, "Sex", covariates)
form <- as.formula(paste("~", paste(rhs, collapse = " + "), "+ (1|participant_id)"))

vobj <- voomWithDreamWeights(dge, form, data)
fit <- eBayes(dream(vobj, form, data))

res <- topTable(fit, coef = "SexFemale", number = Inf, sort.by = "P") %>%
  tibble::rownames_to_column("gene")

write.table(
  res,
  gzfile(file.path(outdir, "DGE_sex.dream.sva.txt.gz")),
  quote = FALSE, sep = "\t", row.names = FALSE
)
