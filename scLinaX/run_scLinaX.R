#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(scLinaX)
  library(readr)
  library(dplyr)
})

## =========================================
## Script: run_scLinaX_summary.R
##
## Usage:
##   Rscript run_scLinaX_summary.R <cohort> [institution]
##
##   cohort:
##     - "MIT_ROSMAP"
##     - "ROSMAP"
##     - "SEA_AD"
##
##   institution (optional):
##     - MIT_ROSMAP: default "R12"
##     - ROSMAP    : default "BroadNYGC_comprehensive"
##     - SEA_AD    : ignored
##
## Environment:
##   XCI_BASE_DIR must point to the base directory that contains
##     - <cohort>/results_*/sclinax (MIT_ROSMAP / ROSMAP)
##     - SEA_AD/results/sclinax
## =========================================

## -----------------------------
## Helper: choose paths by cohort
## -----------------------------

get_base_dir <- function(cohort, institution, xci_base_dir) {
  if (cohort %in% c("MIT_ROSMAP", "ROSMAP")) {
    if (is.na(institution) || institution == "") {
      if (cohort == "MIT_ROSMAP") institution <- "R12"
      if (cohort == "ROSMAP")     institution <- "BroadNYGC_comprehensive"
    }
    base_dir <- file.path(xci_base_dir, cohort, paste0("results_", institution), "sclinax")
  } else if (cohort == "SEA_AD") {
    base_dir <- file.path(xci_base_dir, cohort, "results", "sclinax")
  } else {
    stop("Unknown cohort: ", cohort)
  }
  base_dir
}

get_annotation_path <- function(cohort, xci_base_dir) {
  if (cohort == "MIT_ROSMAP") {
    file.path(xci_base_dir, "MIT_ROSMAP", "annotations",
              "annotations.individualID.remove_duplicates.MIT_ROSMAP.txt.gz")
  } else if (cohort == "ROSMAP") {
    file.path(xci_base_dir, "ROSMAP", "annotations",
              "annotations.ROSMAP.txt.gz")
  } else if (cohort == "SEA_AD") {
    file.path(xci_base_dir, "SEA_AD", "annotations",
              "annotations.SEA_AD.txt.gz")
  } else {
    stop("Unknown cohort: ", cohort)
  }
}

get_notin_sample_table_path <- function(cohort, xci_base_dir) {
  if (cohort == "ROSMAP") {
    file.path(xci_base_dir, "ROSMAP", "data",
              "sample_table_file.female.Broad.notin_MIT_ROSMAP.tsv")
  } else {
    NA_character_
  }
}

get_sea_ad_manifest_path <- function(xci_base_dir) {
  file.path(xci_base_dir, "SEA_AD", "data", "SYNAPSE_METADATA_MANIFEST.tsv")
}

## -----------------------------
## Helper: build only-annotated ASE file if needed
## -----------------------------

prepare_multiome_ase <- function(cohort, only_annotated, base_dir, multiome_annotation) {
  setwd(base_dir)

  if (only_annotated) {
    multiome_ase_path <- "QC_passed_SNP_df.only_annotated.tsv.gz"
    if (!file.exists(multiome_ase_path)) {
      raw_path <- "QC_passed_SNP_df.tsv.gz"
      if (!file.exists(raw_path)) {
        stop("Raw ASE file not found: ", raw_path)
      }
      ase_df <- read_tsv(raw_path, show_col_types = FALSE)

      if (cohort == "MIT_ROSMAP") {
        ase_df$cell_barcode <- paste(ase_df$Sample_ID, ase_df$cell_barcode, sep = "_")
      }

      ase_df <- ase_df %>%
        filter(cell_barcode %in% multiome_annotation$cell_barcode)

      write.table(
        ase_df,
        file = gzfile(multiome_ase_path),
        sep = "\t", quote = FALSE, row.names = FALSE
      )
    }
  } else {
    multiome_ase_path <- "QC_passed_SNP_df.tsv.gz"
  }

  if (!file.exists(multiome_ase_path)) {
    stop("ASE file not found: ", multiome_ase_path)
  }

  read_tsv(multiome_ase_path, show_col_types = FALSE)
}

## -----------------------------
## Main
## -----------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript run_scLinaX_summary.R <cohort> [institution]")
}

cohort      <- as.character(args[1])
institution <- if (length(args) >= 2) as.character(args[2]) else NA_character_

if (!cohort %in% c("MIT_ROSMAP", "ROSMAP", "SEA_AD")) {
  stop("cohort must be one of: MIT_ROSMAP, ROSMAP, SEA_AD")
}

xci_base_dir <- Sys.getenv("XCI_BASE_DIR")
if (xci_base_dir == "") {
  stop("Please set XCI_BASE_DIR environment variable.")
}

base_dir           <- get_base_dir(cohort, institution, xci_base_dir)
annotation_path    <- get_annotation_path(cohort, xci_base_dir)
notin_sample_table <- get_notin_sample_table_path(cohort, xci_base_dir)

message("Cohort: ", cohort)
message("Base dir: ", base_dir)
message("Annotation: ", annotation_path)

if (!dir.exists(base_dir)) {
  stop("Base directory does not exist: ", base_dir)
}

setwd(base_dir)

## Flags per cohort (can be edited if needed)
only_annotated     <- TRUE
notin_MIT_ROSMAP   <- (cohort == "ROSMAP")
separate_replicates <- FALSE
if (cohort == "SEA_AD") {
  separate_replicates <- FALSE
  notin_MIT_ROSMAP   <- FALSE
}

data("XCI_ref")
data("AIDA_QCREF")

multiome_annotation <- read_tsv(annotation_path, show_col_types = FALSE)

multiome_ASE_df <- prepare_multiome_ase(
  cohort              = cohort,
  only_annotated      = only_annotated,
  base_dir            = base_dir,
  multiome_annotation = multiome_annotation
)

if (notin_MIT_ROSMAP) {
  if (is.na(notin_sample_table)) {
    stop("notin_MIT_ROSMAP is TRUE but no sample_table path is defined for this cohort.")
  }
  message("Filtering samples not in MIT_ROSMAP using: ", notin_sample_table)
  sample_table <- read_tsv(notin_sample_table, show_col_types = FALSE)
  multiome_ASE_df <- multiome_ASE_df %>%
    filter(Sample_ID %in% sample_table$sample_name)
}

if (cohort == "SEA_AD") {
  if (!separate_replicates) {
    manifest_path <- get_sea_ad_manifest_path(xci_base_dir)
    message("SEA_AD manifest: ", manifest_path)

    metadata_snRNAseq <- read.delim(
      manifest_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE
    )

    metadata_snRNAseq <- metadata_snRNAseq %>%
      distinct(specimenID, individualID, .keep_all = TRUE)

    multiome_ASE_df <- multiome_ASE_df %>%
      left_join(
        metadata_snRNAseq %>% select(specimenID, individualID),
        by = c("Sample_ID" = "specimenID")
      ) %>%
      mutate(Sample_ID = coalesce(individualID, Sample_ID)) %>%
      select(-individualID)
  }
}

if (cohort %in% c("MIT_ROSMAP", "ROSMAP")) {
  if (is.na(institution) || institution == "") {
    if (cohort == "MIT_ROSMAP") institution <- "R12"
    if (cohort == "ROSMAP")     institution <- "BroadNYGC_comprehensive"
  }
  output_prefix <- paste(cohort, institution, sep = "_")
} else {
  output_prefix <- cohort
}

if (only_annotated) {
  output_prefix <- paste(output_prefix, "only_annotated", sep = ".")
}
if (notin_MIT_ROSMAP) {
  output_prefix <- paste(output_prefix, "notin_MIT_ROSMAP", sep = ".")
}
if (cohort == "SEA_AD" && separate_replicates) {
  output_prefix <- paste(output_prefix, "separate_replicates", sep = ".")
}

message("Output prefix: ", output_prefix)

scLinaX_res <- run_scLinaX(
  ASE_df                  = multiome_ASE_df,
  XCI_ref                 = XCI_ref,
  QCREF                   = AIDA_QCREF,
  Inactive_Gene_ratio_THR = 0.05,
  SNP_DETECTION_DP        = 30,
  SNP_DETECTION_MAF       = 0.1,
  QC_total_allele_THR     = 10,
  HE_allele_cell_number_THR = 50,
  REMOVE_ESCAPE           = TRUE,
  PVAL_THR                = 0.01,
  RHO_THR                 = 0.5
)

saveRDS(scLinaX_res, file = paste0(output_prefix, ".sclinax_results.rds"))

PBMC_summary <- summarize_scLinaX(
  scLinaX_res,
  QC_total_allele_THR = 10,
  Annotation          = NULL
)
write.table(
  PBMC_summary,
  paste0(output_prefix, ".PBMC_summary.txt"),
  sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE
)

cell_type_summary <- summarize_scLinaX(
  scLinaX_res,
  QC_total_allele_THR = 10,
  Annotation          = multiome_annotation
)
write.table(
  cell_type_summary,
  paste0(output_prefix, ".cell_type_summary.txt"),
  sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE
)

message("Done. Generated:")
message("  ", paste0(output_prefix, ".sclinax_results.rds"))
message("  ", paste0(output_prefix, ".PBMC_summary.txt"))
message("  ", paste0(output_prefix, ".cell_type_summary.txt"))
