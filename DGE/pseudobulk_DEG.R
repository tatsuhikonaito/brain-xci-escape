#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(DESeq2)
  library(sva)
})

## =========================================
## Script: pseudobulk_DEG.R
##
## Usage:
##   Rscript pseudobulk_DEG.R \
##       <cohort> <celltype> <analysis> <pheno_definition> [adjust_batch]
##
##   cohort:
##     - "MIT_ROSMAP"
##     - "ROSMAP"
##     - "SEA_AD"
##
##   celltype:
##     - e.g. "MG", "Ast", "all"
##
##   analysis:
##     - "pheno"  : run pheno-based DEG only
##     - "sex"    : run sex effect (fsex) only
##
##   pheno_definition (used only when analysis includes "pheno"):
##     - MIT_ROSMAP / ROSMAP:
##         "AD_else", "AD_NCI", "ADMCI_NCI"
##     - SEA_AD:
##         "AD_else", "AD_NCI"
##
##   adjust_batch (optional):
##     - TRUE/FALSE (default TRUE)
##
## Required environment variables:
##   XCI_BASE_DIR      : base dir for XCI project (cohort subdirs inside)
##   BIGBRAIN_BASE_DIR : base dir for bigbrain data (cohort subdirs inside)
##
## Optional environment variables:
##   SEA_AD_METADATA_PATH  : path to SEA-AD_individual_metadata.csv
##   ROSMAP_CLINICAL_PATH  : path to ROSMAP_clinical.csv
##
## Output examples:
##   <dir_base>/DEG/results/
##     <cohort>_<celltype>.<test>.DEG.<pheno_def>.<sex_group>.pseudobulk.DESeq2[.sva].txt.gz
##     <cohort>_<celltype>.<test>.DEG_sex.pseudobulk.DESeq2[.sva].txt.gz
## =========================================

## -----------------------------
## Helper: pheno-based DEG
## -----------------------------

run_pseudobulk_pheno_deseq <- function(genes_counts,
                                       clinical_filtered,
                                       covariates_base,
                                       tests,
                                       cohort,
                                       celltype,
                                       pheno_definition,
                                       dir_base,
                                       adjust_batch = TRUE) {
  for (test in tests) {
    for (sex_group in c("male", "female", "both_sexes")) {
      covariates <- covariates_base
      message("Pheno-DE | test=", test, " | sex_group=", sex_group)

      if (test == "all") {
        include <- !is.na(clinical_filtered$fsex) &
          apply(!is.na(clinical_filtered[, c(covariates, "pheno"), drop = FALSE]), 1, all)
      } else if (test == "notin_MIT_ROSMAP") {
        samples_notin_other <- scan(
          file.path("data", "sample_list.ROSMAP.notin_MIT_ROSMAP.tsv"),
          what = character(), quiet = TRUE
        )
        include <- !is.na(clinical_filtered$fsex) &
          apply(!is.na(clinical_filtered[, c(covariates, "pheno"), drop = FALSE]), 1, all) &
          clinical_filtered$individualID %in% samples_notin_other
      } else {
        stop("Unknown test: ", test)
      }

      if (sex_group == "male") {
        include <- include & (clinical_filtered$fsex == 0)
      } else if (sex_group == "female") {
        include <- include & (clinical_filtered$fsex == 1)
      } else if (sex_group == "both_sexes") {
        covariates <- c(covariates, "fsex")
      } else {
        stop("Unknown sex_group: ", sex_group)
      }

      design_formula <- as.formula(
        paste("~ pheno +", paste(covariates, collapse = " + "))
      )

      dds <- DESeqDataSetFromMatrix(
        countData = genes_counts[, include, drop = FALSE],
        colData   = clinical_filtered[include, ],
        design    = design_formula
      )

      dds <- dds[rowSums(counts(dds)) > 10, ]

      if (adjust_batch) {
        vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
        dat <- assay(vsd)
        mod  <- model.matrix(design_formula, data = colData(dds))
        mod0 <- model.matrix(~ 1, data = colData(dds))
        nsv_est <- num.sv(dat, mod, method = "leek")

        if (nsv_est > 0) {
          message("Pheno-DE | SVA n.sv = ", nsv_est)
          svobj <- svaseq(dat, mod, mod0, n.sv = nsv_est)
          if (!is.null(svobj$sv)) {
            num_sv <- ncol(svobj$sv)
            for (i in seq_len(num_sv)) {
              colData(dds)[, paste0("SV", i)] <- svobj$sv[, i]
            }
            sv_terms   <- paste(paste0("SV", seq_len(num_sv)), collapse = " + ")
            new_design <- as.formula(
              paste("~", sv_terms, "+ pheno +", paste(covariates, collapse = " + "))
            )
            design(dds) <- new_design
          } else {
            warning("Pheno-DE | SVA returned NULL surrogate variables.")
          }
        } else {
          message("Pheno-DE | n.sv = 0, using original design.")
        }
      }

      dds <- DESeq(dds)

      res_pheno    <- results(dds, name = "pheno", cooksCutoff = FALSE)
      res_pheno_df <- as.data.frame(res_pheno)
      colnames(res_pheno_df) <- paste0(colnames(res_pheno_df), "_pheno")

      if (sex_group == "both_sexes") {
        res_fsex    <- results(dds, name = "fsex", cooksCutoff = FALSE)
        res_fsex_df <- as.data.frame(res_fsex)
        colnames(res_fsex_df) <- paste0(colnames(res_fsex_df), "_fsex")
        res_merge <- merge(res_fsex_df, res_pheno_df, by = "row.names", all = TRUE)
        rownames(res_merge) <- res_merge$Row.names
        res_merge$Row.names <- NULL
      } else {
        res_merge <- res_pheno_df
        rownames(res_merge) <- rownames(res_pheno_df)
      }

      if (adjust_batch) {
        outfile <- file.path(
          dir_base, "DEG", "results",
          paste0(
            cohort, "_", celltype, ".", test,
            ".DEG.", pheno_definition, ".", sex_group,
            ".pseudobulk.DESeq2.sva.txt.gz"
          )
        )
      } else {
        outfile <- file.path(
          dir_base, "DEG", "results",
          paste0(
            cohort, "_", celltype, ".", test,
            ".DEG.", pheno_definition, ".", sex_group,
            ".pseudobulk.DESeq2.txt.gz"
          )
        )
      }

      dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
      write.table(res_merge, gzfile(outfile), quote = FALSE, sep = "\t")
      message("Pheno-DE | Saved: ", outfile)
    }
  }
}

## -----------------------------
## Helper: sex-effect DEG
## -----------------------------

run_pseudobulk_sex_deseq <- function(genes_counts,
                                     clinical_filtered,
                                     covariates_base,
                                     tests,
                                     cohort,
                                     celltype,
                                     dir_base,
                                     adjust_batch = TRUE) {
  for (test in tests) {
    covariates <- covariates_base
    message("Sex-DE | test=", test)

    if (test == "all") {
      include <- !is.na(clinical_filtered$fsex) &
        apply(!is.na(clinical_filtered[, covariates, drop = FALSE]), 1, all)
    } else if (test == "notin_MIT_ROSMAP") {
      samples_notin_other <- scan(
        file.path("data", "sample_list.ROSMAP.notin_MIT_ROSMAP.tsv"),
        what = character(), quiet = TRUE
      )
      include <- !is.na(clinical_filtered$fsex) &
        apply(!is.na(clinical_filtered[, covariates, drop = FALSE]), 1, all) &
        clinical_filtered$individualID %in% samples_notin_other
    } else {
      stop("Unknown test: ", test)
    }

    design_formula <- as.formula(
      paste("~ fsex +", paste(covariates, collapse = " + "))
    )

    dds <- DESeqDataSetFromMatrix(
      countData = genes_counts[, include, drop = FALSE],
      colData   = clinical_filtered[include, ],
      design    = design_formula
    )

    dds <- dds[rowSums(counts(dds)) > 10, ]

    if (adjust_batch) {
      vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
      dat <- assay(vsd)
      mod  <- model.matrix(design_formula, data = colData(dds))
      mod0 <- model.matrix(~ 1, data = colData(dds))
      nsv_est <- num.sv(dat, mod, method = "leek")

      if (nsv_est > 0) {
        message("Sex-DE | SVA n.sv = ", nsv_est)
        svobj <- svaseq(dat, mod, mod0, n.sv = nsv_est)
        if (!is.null(svobj$sv)) {
          num_sv <- ncol(svobj$sv)
          for (i in seq_len(num_sv)) {
            colData(dds)[, paste0("SV", i)] <- svobj$sv[, i]
          }
          sv_terms   <- paste(paste0("SV", seq_len(num_sv)), collapse = " + ")
          new_design <- as.formula(
            paste("~", sv_terms, "+ fsex +", paste(covariates, collapse = " + "))
          )
          design(dds) <- new_design
        } else {
          warning("Sex-DE | SVA returned NULL surrogate variables.")
        }
      } else {
        message("Sex-DE | n.sv = 0, using original design.")
      }
    }

    dds <- DESeq(dds)

    res_fsex    <- results(dds, name = "fsex", cooksCutoff = FALSE)
    res_fsex_df <- as.data.frame(res_fsex)

    if (adjust_batch) {
      outfile <- file.path(
        dir_base, "DEG", "results",
        paste0(
          cohort, "_", celltype, ".", test,
          ".DEG_sex.pseudobulk.DESeq2.sva.txt.gz"
        )
      )
    } else {
      outfile <- file.path(
        dir_base, "DEG", "results",
        paste0(
          cohort, "_", celltype, ".", test,
          ".DEG_sex.pseudobulk.DESeq2.txt.gz"
        )
      )
    }

    dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
    write.table(res_fsex_df, gzfile(outfile), quote = FALSE, sep = "\t")
    message("Sex-DE | Saved: ", outfile)
  }
}

## -----------------------------
## Main
## -----------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript pseudobulk_DEG.R <cohort> <celltype> <analysis> <pheno_definition> [adjust_batch]")
}

cohort           <- as.character(args[1])
celltype         <- as.character(args[2])
analysis         <- as.character(args[3])  # "pheno", "sex"
pheno_definition <- as.character(args[4])
adjust_batch     <- if (length(args) >= 5) as.logical(args[5]) else TRUE

if (!cohort %in% c("MIT_ROSMAP", "ROSMAP", "SEA_AD")) {
  stop("Unknown cohort: ", cohort)
}
if (!analysis %in% c("pheno", "sex")) {
  stop("analysis must be one of: pheno or sex")
}

xci_base      <- Sys.getenv("XCI_BASE_DIR")
bigbrain_base <- Sys.getenv("BIGBRAIN_BASE_DIR")
if (xci_base == "" || bigbrain_base == "") {
  stop("Please set XCI_BASE_DIR and BIGBRAIN_BASE_DIR environment variables.")
}

dir_base <- file.path(xci_base, cohort)
dir_data <- file.path(bigbrain_base, cohort)
setwd(dir_base)

tests <- if (cohort == "ROSMAP") c("all", "notin_MIT_ROSMAP") else "all"

message("Cohort: ", cohort)
message("Cell type: ", celltype)
message("Analysis: ", analysis)
message("Pheno definition: ", pheno_definition, " (used only if analysis includes pheno)")
message("Adjust batch (SVA): ", adjust_batch)
message("XCI base dir: ", dir_base)
message("Data dir: ", dir_data)

## =============================
## Clinical data loading
## =============================

if (cohort == "MIT_ROSMAP" || cohort == "ROSMAP") {
  rosmap_clinical_path <- Sys.getenv("ROSMAP_CLINICAL_PATH")
  if (rosmap_clinical_path == "") {
    rosmap_clinical_path <- file.path(
      bigbrain_base, "ROSMAP", "metadata", "ROSMAP_clinical.csv"
    )
  }
  message("Clinical file (ROSMAP-based): ", rosmap_clinical_path)
  clinical <- read_csv(rosmap_clinical_path, show_col_types = FALSE)

  clinical$cogdx <- as.character(clinical$cogdx)
  clinical <- clinical %>%
    mutate(
      fsex = case_when(
        !is.na(msex) & msex == 1 ~ 0,
        !is.na(msex) & msex == 0 ~ 1,
        TRUE ~ NA_real_
      ),
      apoe4 = case_when(
        is.na(apoe_genotype) ~ NA_real_,
        grepl("4", as.character(apoe_genotype)) ~ 1,
        TRUE ~ 0
      )
    )
}

if (cohort == "SEA_AD") {
  sea_ad_meta_path <- Sys.getenv("SEA_AD_METADATA_PATH")
  if (sea_ad_meta_path == "") {
    sea_ad_meta_path <- file.path(
      xci_base, "SEA_AD", "data", "SEA-AD_individual_metadata.csv"
    )
  }
  message("Clinical file (SEA_AD): ", sea_ad_meta_path)
  clinical <- read_csv(sea_ad_meta_path, show_col_types = FALSE)

  clinical <- clinical %>%
    mutate(
      apoe4 = case_when(
        !is.na(apoe4Status) & apoe4Status == "Y" ~ 1,
        !is.na(apoe4Status) & apoe4Status == "N" ~ 0,
        TRUE ~ NA_real_
      )
    )
}

## =============================
## Pseudobulk counts
## =============================

message("Loading pseudobulk data...")
if (celltype == "all") {
  load(file.path("pseudobulk", paste0(cohort, "_all_RNA_pseudobulk_EUR.RData")))
} else {
  load(file.path(
    dir_data, "analysis", "snRNAseq", "seurat",
    paste0(cohort, "_", celltype, "_RNA_psuedobulk_EUR.RData")
  ))
}

gene_names   <- as.vector(genes_counts[[1]])
genes_counts <- as.matrix(genes_counts[, -1, drop = FALSE])
rownames(genes_counts) <- gene_names

## =============================
## Align clinical and counts
## =============================

clinical_filtered <- clinical[clinical$individualID %in% colnames(genes_counts), ]
clinical_filtered <- clinical_filtered[
  match(colnames(genes_counts), clinical_filtered$individualID),
]

## =============================
## Cohort-specific formatting
## =============================

if (cohort == "MIT_ROSMAP" || cohort == "ROSMAP") {
  clinical_filtered$cogdx     <- as.character(clinical_filtered$cogdx)
  clinical_filtered$Study     <- as.character(clinical_filtered$Study)
  clinical_filtered$age_death <- as.numeric(
    gsub("90\\+", "90", clinical_filtered$age_death)
  )

  ## pheno coding (used only if analysis includes pheno)
  if (analysis == "pheno") {
    if (pheno_definition == "AD_else") {
      clinical_filtered$pheno <- ifelse(
        clinical_filtered$cogdx %in% c("4", "5"), 1,
        ifelse(clinical_filtered$cogdx %in% c("1", "2", "3", "6"), 0, NA)
      )
    } else if (pheno_definition == "AD_NCI") {
      clinical_filtered$pheno <- ifelse(
        clinical_filtered$cogdx %in% c("4", "5"), 1,
        ifelse(clinical_filtered$cogdx %in% "1", 0, NA)
      )
    } else if (pheno_definition == "ADMCI_NCI") {
      clinical_filtered$pheno <- ifelse(
        clinical_filtered$cogdx %in% c("2", "3", "4", "5"), 1,
        ifelse(clinical_filtered$cogdx %in% "1", 0, NA)
      )
    } else {
      stop("Unknown pheno definition for MIT_ROSMAP/ROSMAP: ", pheno_definition)
    }
  }

  covariates_base <- c("age_death", "pmi", "apoe4", "Study")
  if (analysis == "sex") {
    covariates_base <- c(covariates_base, "cogdx")
  }
}

if (cohort == "SEA_AD") {
  clinical_filtered$ageDeath <- as.numeric(
    gsub("90\\+", "90", clinical_filtered$ageDeath)
  )

  clinical_filtered <- clinical_filtered %>%
    mutate(
      pheno = case_when(
        diagnosis %in% c(
          "\"Alzheimer disease,dementia\"",
          "\"Alzheimer disease,dementia,Lewy body disease\"",
          "\"Alzheimer disease,dementia,vascular dementia\""
        ) ~ "AD",
        diagnosis == "\"control,no cognitive impairment\"" ~ "control",
        diagnosis %in% c(
          "\"mild cognitive impairment\"",
          "\"dementia\"",
          "\"dementia,Brain Cancer\"",
          "\"dementia,multiple system atrophy\"",
          "\"dementia,Parkinson's disease\"",
          "\"dementia,vascular dementia\""
        ) ~ "others",
        diagnosis == "Not Applicable" ~ "unknown",
        TRUE ~ NA_character_
      ),
      fsex = case_when(
        sex == "male" ~ 0,
        sex == "female" ~ 1,
        TRUE ~ NA_real_
      )
    )

  if (analysis == "pheno") {
    if (pheno_definition == "AD_else") {
      clinical_filtered$pheno <- ifelse(
        clinical_filtered$pheno %in% "AD", 1,
        ifelse(clinical_filtered$pheno %in% c("control", "others", "unknown"), 0, NA)
      )
    } else if (pheno_definition == "AD_NCI") {
      clinical_filtered$pheno <- ifelse(
        clinical_filtered$pheno %in% "AD", 1,
        ifelse(clinical_filtered$pheno %in% "control", 0, NA)
      )
    } else {
      stop("Unknown pheno definition for SEA_AD: ", pheno_definition)
    }
  }

  covariates_base <- c("ageDeath", "pmi", "apoe4")
  if (analysis == "sex") {
    covariates_base <- c(covariates_base, "pheno")
  }
}

## For MIT_ROSMAP / ROSMAP, fsex already defined in clinical; align filtered if needed
if (cohort == "MIT_ROSMAP" || cohort == "ROSMAP") {
  if (!"fsex" %in% colnames(clinical_filtered)) {
    stop("fsex is missing in clinical_filtered for cohort: ", cohort)
  }
}

## =============================
## Run analyses
## =============================

if (analysis == "pheno") {
  run_pseudobulk_pheno_deseq(
    genes_counts      = genes_counts,
    clinical_filtered = clinical_filtered,
    covariates_base   = covariates_base,
    tests             = tests,
    cohort            = cohort,
    celltype          = celltype,
    pheno_definition  = pheno_definition,
    dir_base          = dir_base,
    adjust_batch      = adjust_batch
  )
}

else if (analysis == "sex") {
  run_pseudobulk_sex_deseq(
    genes_counts      = genes_counts,
    clinical_filtered = clinical_filtered,
    covariates_base   = covariates_base,
    tests             = tests,
    cohort            = cohort,
    celltype          = celltype,
    dir_base          = dir_base,
    adjust_batch      = adjust_batch
  )
}
