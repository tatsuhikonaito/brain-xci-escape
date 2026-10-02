#!/usr/bin/env bash
# Single-sample variant calling and two-tier filtering.
# Run from the repository root. The sites VCF must match this sample's ancestry.
set -euo pipefail

if [[ $# -ne 4 ]]; then
  echo "Usage: bash bulk_microglia/02_prepare_XCIR_input.sh <sample.bam> <reference.fa> <population.sites.vcf.gz> <output_prefix>" >&2
  exit 1
fi
bam="$1"
ref_fa="$2"
sites_vcf="$3"
out_prefix="$4"
[[ -s "$bam" ]] || { echo "BAM not found: $bam" >&2; exit 1; }
mkdir -p "$(dirname "$out_prefix")"

# BAM, reference and sites VCF must use the same contig names.
contigs="$(samtools idxstats "$bam" | cut -f1)"
if grep -x 'chrX' <<< "$contigs" > /dev/null; then
  chr_x="chrX"
elif grep -x 'X' <<< "$contigs" > /dev/null; then
  chr_x="X"
else
  echo "chrX/X not found in BAM: $bam" >&2
  exit 1
fi
# Prepare the shared sites index before launching concurrent samples.
if [[ ! -f "${sites_vcf}.tbi" ]]; then
  tabix -p vcf "$sites_vcf"
fi

MIN_MAPQ=20
MIN_BQ=20
MAX_DP=100000
T1_MIN_DP=6
T1_MIN_AD_EACH=1
T2_MIN_DP=20
T2_MIN_AD_EACH=3

prefix="${out_prefix}.${chr_x}"
all_raw="${prefix}.all.raw.vcf.gz"
tier1_raw="${prefix}.tier1.raw.vcf.gz"
tier2_raw="${prefix}.tier2.raw.vcf.gz"
tier1_filt="${prefix}.tier1.filt.vcf.gz"
tier2_filt="${prefix}.tier2.filt.vcf.gz"
combined_vcf="${prefix}.XCIR_input.vcf.gz"

# MAX_DP is applied to the VCF, not to the mpileup depth setting.
bcftools mpileup -f "$ref_fa" -r "$chr_x" \
  -q "$MIN_MAPQ" -Q "$MIN_BQ" -a FORMAT/AD,FORMAT/DP -Ou "$bam" \
| bcftools call -mv -Ou \
| bcftools norm -f "$ref_fa" -Oz -o "$all_raw"
tabix -p vcf "$all_raw"

bcftools isec -n=2 -w1 "$all_raw" "$sites_vcf" -Oz -o "$tier1_raw"
bcftools isec -n=1 -w1 "$all_raw" "$sites_vcf" -Oz -o "$tier2_raw"
tabix -p vcf "$tier1_raw"
tabix -p vcf "$tier2_raw"

bcftools view -m2 -M2 -v snps "$tier1_raw" -Ou \
| bcftools filter -e "FMT/DP<${T1_MIN_DP} || FMT/DP>${MAX_DP}" -Ou \
| bcftools filter -e "FMT/AD[0:0]<${T1_MIN_AD_EACH} || FMT/AD[0:1]<${T1_MIN_AD_EACH}" \
  -Oz -o "$tier1_filt"
tabix -p vcf "$tier1_filt"

bcftools view -m2 -M2 -v snps "$tier2_raw" -Ou \
| bcftools filter -e "FMT/DP<${T2_MIN_DP} || FMT/DP>${MAX_DP}" -Ou \
| bcftools filter -e "FMT/AD[0:0]<${T2_MIN_AD_EACH} || FMT/AD[0:1]<${T2_MIN_AD_EACH}" \
  -Oz -o "$tier2_filt"
tabix -p vcf "$tier2_filt"

bcftools concat -a "$tier1_filt" "$tier2_filt" -Oz -o "$combined_vcf"
tabix -p vcf "$combined_vcf"

echo "XCIR input: $combined_vcf"
echo "Tier1 SNPs: $(bcftools view -H "$tier1_filt" | wc -l)"
echo "Tier2 SNPs: $(bcftools view -H "$tier2_filt" | wc -l)"
