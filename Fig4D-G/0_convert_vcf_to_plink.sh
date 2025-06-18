#!/bin/bash

#===============================================================================
# Author: Marco A. Coelho @ Heitman Lab, Duke University
# Date: 2025-05-28
# Version: 1.0
#
# Description:
#   Converts a multi-sample VCF from Snippy into PLINK binary format
#   (.bed, .bim, .fam) for PCA, LD decay, and other population genomics analyses.
#
# Requirements:
#   - Conda environment: popgen_env
#   - Tool: PLINK v1.9 or newer (https://www.cog-genomics.org/plink/)
#
# Usage:
#   conda activate popgen_env
#   bash 0_convert_vcf_to_plink.sh
#
# Input:
#   core.vcf            # Multi-sample VCF file
#
# Output:
#   core.bed            # PLINK binary genotype table
#   core.bim            # SNP information
#   core.fam            # Sample information
#===============================================================================

set -euo pipefail

# --- INPUT ---
VCF_FILE="core.vcf"
OUT_PREFIX="core"

# --- CHECK INPUT ---
if [[ ! -f "$VCF_FILE" ]]; then
    echo "Error: Input VCF not found at $VCF_FILE"
    exit 1
fi

# --- CHECK PLINK ---
if ! command -v plink &> /dev/null; then
    echo "Error: PLINK not found in PATH. Please activate popgen_env."
    exit 1
fi

# --- CONVERT ---
echo "Converting $VCF_FILE to PLINK binary format..."
plink --vcf "$VCF_FILE" --make-bed --allow-extra-chr --out "$OUT_PREFIX"

# --- VERIFY OUTPUT ---
if [[ -f "${OUT_PREFIX}.bed" && -f "${OUT_PREFIX}.bim" && -f "${OUT_PREFIX}.fam" ]]; then
    echo "PLINK conversion successful: ${OUT_PREFIX}.bed/.bim/.fam"
else
    echo "Conversion failed: Expected PLINK files not found."
    exit 1
fi
