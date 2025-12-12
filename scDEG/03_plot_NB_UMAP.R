#!/usr/bin/env Rscript

# Plot NB-based sex effect (Beta_total) and gene expression on UMAP

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrastr)
  library(viridis)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7) {
  stop("Usage: Rscript plot_NB_UMAP.R COV RES CELL_TYPE_INFO UMAP UMI_Count GENE OUTNAME")
}

COV       <- args[1]
RES       <- args[2]
UMAP_FILE <- args[4]
UMI_FILE  <- args[5]
GENE      <- args[6]
OUTNAME   <- args[7]

## 1. metadata + NB results → Beta_total
cov <- readr::read_tsv(COV, show_col_types = FALSE)
res <- readr::read_tsv(RES, show_col_types = FALSE) |>
  dplyr::filter(Gene == GENE, stringr::str_detect(Variable, "fsex"))
if (nrow(res) == 0) stop("No NB result found for gene: ", GENE)

beta_v <- res$BETA
names(beta_v) <- res$Variable

hmPC_cols  <- paste0("hmPC_", 1:10)
beta_int   <- beta_v[paste0("fsex:hmPC_", 1:10)]

if (!all(hmPC_cols %in% colnames(cov))) stop("hmPC_1–hmPC_10 missing in metadata.")
if (!"fsex" %in% colnames(cov))         stop("fsex column is missing in metadata.")
if (!"barcode" %in% colnames(cov))      stop("barcode column is missing in metadata.")

cov <- cov |>
  mutate(Beta_total = as.numeric(beta_v["fsex"] + as.matrix(across(all_of(hmPC_cols))) %*% beta_int))

beta_umap <- cov |>
  dplyr::select(barcode, Beta_total) |>
  dplyr::filter(!is.na(Beta_total))

## 2. UMAP
umap <- readr::read_tsv(UMAP_FILE, show_col_types = FALSE)
if (ncol(umap) < 3) stop("UMAP file must have ≥3 columns: barcode, UMAP1, UMAP2")
colnames(umap)[1:3] <- c("barcode", "UMAP1", "UMAP2")

beta_umap <- beta_umap |>
  inner_join(umap, by = "barcode")

## 3. Beta_total on UMAP (viridis)
p99     <- stats::quantile(beta_umap$Beta_total, probs = 0.99, na.rm = TRUE)
max_abs <- abs(p99)

p_beta <- ggplot(beta_umap, aes(x = UMAP1, y = UMAP2, colour = Beta_total)) +
  geom_point_rast(size = 0.01, alpha = 0.5) +
  scale_color_viridis(option = "viridis",
                      limits = c(-max_abs, max_abs),
                      name   = expression(beta[total])) +
  theme_test() +
  theme(
    aspect.ratio = 1,
    panel.border = element_blank(),
    axis.line    = element_blank(),
    axis.ticks   = element_blank(),
    axis.text    = element_blank(),
    axis.title   = element_blank(),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 12)
  )

ggsave(paste0(OUTNAME, "_Beta_total_on_UMAP_viridis.pdf"),
       p_beta, width = 4.8, height = 3.6)

## 4. normalized UMI for target gene on UMAP
umi_raw <- readr::read_tsv(UMI_FILE, show_col_types = FALSE)
if (!"barcode" %in% colnames(umi_raw)) stop("UMI file must have column 'barcode'")
if (!GENE %in% colnames(umi_raw))      stop("Gene ", GENE, " not found in UMI file")

umi_expr <- umi_raw |>
  mutate(total_UMI = rowSums(across(-barcode))) |>
  transmute(
    barcode  = barcode,
    EXP      = .data[[GENE]],
    CP10K    = (EXP / total_UMI) * 1e4,
    log_norm = log1p(CP10K)
  )

umi_umap <- umi_expr |>
  inner_join(umap, by = "barcode")

p_expr <- ggplot(umi_umap, aes(x = UMAP1, y = UMAP2, colour = log_norm)) +
  geom_point_rast(size = 0.01, alpha = 0.5) +
  scale_color_viridis(option = "magma",
                      name   = "Normalized\nUMI (log1p(CP10K))") +
  theme_test() +
  theme(
    aspect.ratio = 1,
    panel.border = element_blank(),
    axis.line    = element_blank(),
    axis.ticks   = element_blank(),
    axis.text    = element_blank(),
    axis.title   = element_blank(),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 12)
  )

ggsave(paste0(OUTNAME, "_Norm_UMI_on_UMAP.pdf"),
       p_expr, width = 4.8, height = 3.6)
