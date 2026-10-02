#!/usr/bin/env Rscript

# NB-based sex association test (MIT_ROSMAP)

suppressPackageStartupMessages({
  library(lme4)
  library(Matrix)
  library(tidyverse)
})

# -----------------------------
# Settings
# -----------------------------
cohort   <- "MIT_ROSMAP"
base_dir <- "/path/to/XCI"   # <-- change here

meta_file <- file.path(base_dir, cohort, "data",
                       paste0(cohort, ".metadata.txt.gz"))
expr_file <- file.path(base_dir, cohort, "data",
                       paste0(cohort, ".expression.chrX.txt.gz"))
out_dir   <- file.path(base_dir, cohort, "NB", "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

samples <- readLines(
  file.path(base_dir, cohort, "data",
            paste0("sample_list.", cohort, ".pseudobulk.txt"))
)

# covariates used in the NB model (MIT_ROSMAP)
list_covariates <- c(
  "age", "pmi", "cogdx", "apoe4", "Study", "nUMI",
  stringr::str_c("oriPC_", 1:10)
)

covariates <- c(
  "age", "pmi", "cogdx", "apoe4", "Study", "nUMI",
  stringr::str_c("oriPC_", 1:10),
  "(1 | individualID)"
)

covariates_hmPC <- c(
  covariates,
  stringr::str_c("hmPC_", 1:10)
)

set.seed(9)

# -----------------------------
# Load data
# -----------------------------
meta <- readr::read_tsv(meta_file, show_col_types = FALSE) |>
  dplyr::filter(individualID %in% samples)

expr <- readr::read_tsv(expr_file, show_col_types = FALSE)

gene_name <- colnames(expr)[-1]
# For selected genes: gene_name <- c("LAMP2", "DDX3X")

# -----------------------------
# Result containers
# -----------------------------
SEX_EXP_summary       <- tibble()
SEX_EXP_LRT           <- tibble()
SEX_EXP_PCINT_summary <- tibble()
SEX_EXP_PCINT_LRT     <- tibble()

# -----------------------------
# Per-gene loop
# -----------------------------
for (GENE in gene_name) {
  message("Processing ", GENE)

  test_df <- expr |>
    dplyr::select("barcode", all_of(GENE)) |>
    dplyr::inner_join(meta, by = "barcode") |>
    dplyr::rename(EXP = !!GENE) |>
    dplyr::mutate(cogdx = as.character(cogdx))

  include <- !is.na(test_df$fsex) &
    apply(!is.na(test_df[, list_covariates, drop = FALSE]), 1, all)

  test_df <- test_df[include, , drop = FALSE]

  # ---------------------------
  # 1) Sex main effect
  # ---------------------------
  full_form <- stringr::str_c(c("EXP ~ fsex", covariates), collapse = " + ")
  null_form <- stringr::str_c(c("EXP ~",      covariates), collapse = " + ")

  check <- try({
    full_model <- lme4::glmer.nb(
      formula = as.formula(full_form),
      data    = test_df,
      nAGQ    = 0,
      control = glmerControl(optimizer = "nloptwrap")
    )

    null_model <- lme4::glmer.nb(
      formula = as.formula(null_form),
      data    = test_df,
      nAGQ    = 0,
      control = glmerControl(optimizer = "nloptwrap")
    )

    model_lrt <- as.data.frame(anova(null_model, full_model))
    lrt_df <- tibble(
      Gene  = GENE,
      MODEL = rownames(model_lrt),
      as_tibble(model_lrt)
    )
    SEX_EXP_LRT <- bind_rows(SEX_EXP_LRT, lrt_df)

    coef_df <- as.data.frame(summary(full_model)$coefficients)
    colnames(coef_df) <- c("BETA", "SE", "Z", "P")

    tbl_out <- tibble(
      Gene     = GENE,
      Variable = rownames(coef_df),
      as_tibble(coef_df)
    )
    SEX_EXP_summary <- bind_rows(SEX_EXP_summary, tbl_out)
  })

  if (inherits(check, "try-error")) {
    lrt_df <- tibble(
      Gene        = GENE,
      MODEL       = NA_character_,
      npar        = NA_real_,
      AIC         = NA_real_,
      BIC         = NA_real_,
      logLik      = NA_real_,
      deviance    = NA_real_,
      Chisq       = NA_real_,
      Df          = NA_real_,
      `Pr(>Chisq)`= NA_real_
    )
    SEX_EXP_LRT <- bind_rows(SEX_EXP_LRT, lrt_df)

    tbl_out <- tibble(
      Gene     = GENE,
      Variable = NA_character_,
      BETA     = NA_real_,
      SE       = NA_real_,
      Z        = NA_real_,
      P        = NA_real_
    )
    SEX_EXP_summary <- bind_rows(SEX_EXP_summary, tbl_out)
  }

  # ---------------------------
  # 2) Sex × hmPC interaction
  # ---------------------------
  full_form_int <- stringr::str_c(
    c("EXP ~ fsex", covariates_hmPC,
      stringr::str_c("fsex:hmPC_", 1:10)),
    collapse = " + "
  )
  null_form_int <- stringr::str_c(
    c("EXP ~ fsex", covariates_hmPC),
    collapse = " + "
  )

  check2 <- try({
    full_model_int <- lme4::glmer.nb(
      formula = as.formula(full_form_int),
      data    = test_df,
      nAGQ    = 0,
      control = glmerControl(optimizer = "nloptwrap")
    )

    null_model_int <- lme4::glmer.nb(
      formula = as.formula(null_form_int),
      data    = test_df,
      nAGQ    = 0,
      control = glmerControl(optimizer = "nloptwrap")
    )

    model_lrt2 <- as.data.frame(anova(null_model_int, full_model_int))
    lrt_df2 <- tibble(
      Gene  = GENE,
      MODEL = rownames(model_lrt2),
      as_tibble(model_lrt2)
    )
    SEX_EXP_PCINT_LRT <- bind_rows(SEX_EXP_PCINT_LRT, lrt_df2)

    coef_df2 <- as.data.frame(summary(full_model_int)$coefficients)
    colnames(coef_df2) <- c("BETA", "SE", "Z", "P")

    tbl_out2 <- tibble(
      Gene     = GENE,
      Variable = rownames(coef_df2),
      as_tibble(coef_df2)
    )
    SEX_EXP_PCINT_summary <- bind_rows(SEX_EXP_PCINT_summary, tbl_out2)
  })

  if (inherits(check2, "try-error")) {
    lrt_df2 <- tibble(
      Gene        = GENE,
      MODEL       = NA_character_,
      npar        = NA_real_,
      AIC         = NA_real_,
      BIC         = NA_real_,
      logLik      = NA_real_,
      deviance    = NA_real_,
      Chisq       = NA_real_,
      Df          = NA_real_,
      `Pr(>Chisq)`= NA_real_
    )
    SEX_EXP_PCINT_LRT <- bind_rows(SEX_EXP_PCINT_LRT, lrt_df2)

    tbl_out2 <- tibble(
      Gene     = GENE,
      Variable = NA_character_,
      BETA     = NA_real_,
      SE       = NA_real_,
      Z        = NA_real_,
      P        = NA_real_
    )
    SEX_EXP_PCINT_summary <- bind_rows(SEX_EXP_PCINT_summary, tbl_out2)
  }

  rm(test_df)
}

outprefix <- file.path(out_dir, cohort)

readr::write_tsv(
  SEX_EXP_summary,
  paste0(outprefix, "_Sex_Exp_assoc_NB.summary.tsv.gz")
)
readr::write_tsv(
  SEX_EXP_LRT,
  paste0(outprefix, "_Sex_Exp_assoc_NB.LRT.tsv.gz")
)
readr::write_tsv(
  SEX_EXP_PCINT_summary,
  paste0(outprefix, "_Sex_Exp_PCint_assoc_NB.summary.tsv.gz")
)
readr::write_tsv(
  SEX_EXP_PCINT_LRT,
  paste0(outprefix, "_Sex_Exp_PCint_assoc_NB.LRT.tsv.gz")
)
