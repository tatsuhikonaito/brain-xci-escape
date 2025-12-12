#!/usr/bin/env bash

###############################################################################
# MAGMA gene / gene-set analysis on X chromosome
#
# - Input:
#     * X-only GWAS summary stats (TSV, no header)
#     * MAGMA SNP location file for X (hg38)
#     * MAGMA gene location file for X (hg38)
#     * LD reference panel for X (bfile)
# - Output:
#     * MAGMA gene results
#     * MAGMA gene-set results
#
###############################################################################

########################
# User settings
########################

# Project root
BASE_DIR="/path/to/XCI_project"
cd "${BASE_DIR}"

# Directory for summary statistics
SS_DIR="${BASE_DIR}/SS"

# MAGMA reference directory (bfiles, .loc, etc.)
REF_DIR="/path/to/magma_ref"

# MAGMA binary
MAGMA_BIN="/path/to/magma"

# Gene location file for chrX (hg38, MAGMA format)
GENE_LOC="${REF_DIR}/gencode.v38.chrX.MAGMA.gene.loc"

# X-chromosome reference (hg38, PLINK bfile prefix)
BFILE_X_FEMALE_NONPAR="${REF_DIR}/g1000_eur_X_female_hg38.nonPAR"

# SNP location file for X (hg38, MAGMA snp-loc format)
SNP_LOC="${REF_DIR}/magma_X.snploc.txt"

# Gene set file (MAGMA set annotation format)
GENESET_FILE="${BASE_DIR}/MAGMA/data/magma_genesets.txt"

# Output directory
RESULTS_DIR="${BASE_DIR}/MAGMA/results"
mkdir -p "${RESULTS_DIR}"

########################
# Phenotype settings
########################
PHENO="GCST90444373"

# --- phenotype-specific column mapping (no header TSV) ---
# GCST90444373:
#   8:p_value, 9:n, 10:rs_id
# GCST90449045:
#   8:p_value, 10:rsid, 19:n_total
if [ "${PHENO}" = "GCST90444373" ]; then
  N=1152284
  N_COL=9
else if [ "${PHENO}" = "GCST90449045" ]; then
  N=870171
  N_COL=19
fi


SS_IN="${SS_DIR}/${PHENO}.tsv"
SS_MAGMA="${SS_DIR}/${PHENO}.MAGMA.txt"

###############################################################################
# 1. Prepare MAGMA summary statistics file
#   For MAGMA: SNP (rsid), P, N (total sample size).
###############################################################################

echo -e "SNP\tP\tN" > "${SS_MAGMA}"

awk -F'\t' -v maxN="${N}" -v ncol="${N_COL}" -v pcol="${P_COL}" -v rscol="${RSID_COL}" '
  NR > 1 {
    snp = $rscol
    if (snp == "") next
    if (seen[snp]++) next

    p = $pcol
    n = $(ncol) + 0

    if (n >= 0.95 * maxN) {
      print snp "\t" p "\t" n
    }
  }
' "${SS_IN}" >> "${SS_MAGMA}"


###############################################################################
# 2. MAGMA: annotate SNP → gene
###############################################################################

ANNOT_OUT="${RESULTS_DIR}/magma_X_annot"

"${MAGMA_BIN}" \
  --annotate \
  --snp-loc "${SNP_LOC}" \
  --gene-loc "${GENE_LOC}" \
  --out "${ANNOT_OUT}"

###############################################################################
# 3. MAGMA: gene-based test on X
###############################################################################

GENE_OUT_NONPAR="${RESULTS_DIR}/${PHENO}.magma_X_gene.nonPAR"

"${MAGMA_BIN}" \
  --bfile "${BFILE_X_FEMALE_NONPAR}" \
  --pval "${SS_MAGMA}" N=${N} \
  --gene-annot "${ANNOT_OUT}.genes.annot" \
  --out "${GENE_OUT_NONPAR}"

###############################################################################
# 4. MAGMA: gene-set analysis
###############################################################################

GENESET_OUT_NONPAR="${RESULTS_DIR}/${PHENO}.magma_X_geneset.nonPAR"

"${MAGMA_BIN}" \
  --gene-results "${GENE_OUT_NONPAR}.genes.raw" \
  --set-annot "${GENESET_FILE}" \
  --out "${GENESET_OUT_NONPAR}"

