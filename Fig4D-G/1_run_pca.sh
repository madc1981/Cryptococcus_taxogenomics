#!/bin/bash

#===============================================================================
# Author: Marco A. Coelho @ Heitman Lab, Duke University
# Date: 2025-05-28
# Version: 1.0
#
# Description:
#   Runs Principal Component Analysis (PCA) on PLINK binary genotype files
#   (.bed, .bim, .fam), outputting eigenvectors and eigenvalues for visualization
#   and population structure analysis.
#
# Requirements:
#   - Conda environment: popgen_env
#   - Tool: PLINK v1.9 or newer
#
# Usage:
#   conda activate popgen_env
#   bash 1_run_pca.sh
#
# Input:
#   core.bed, core.bim, core.fam     # PLINK binary genotype files
#
# Output:
#   core_pca.eigenvec                # Principal component coordinates
#   core_pca.eigenval                # Variance explained by each PC
#===============================================================================

set -euo pipefail

# --- INPUT PREFIX ---
INPUT_PREFIX="core"
OUT_PREFIX="core_pca"

# --- CHECK INPUT ---
if [[ ! -f "${INPUT_PREFIX}.bed" ]]; then
    echo "Error: PLINK .bed file not found for input prefix '$INPUT_PREFIX'"
    exit 1
fi

# --- RUN PCA ---
echo "Running PCA with PLINK..."
plink --bfile "$INPUT_PREFIX" --pca --allow-extra-chr --out "$OUT_PREFIX"

# --- CHECK OUTPUT ---
if [[ -f "${OUT_PREFIX}.eigenvec" && -f "${OUT_PREFIX}.eigenval" ]]; then
    echo "PCA complete:"
    echo "  - Eigenvectors: ${OUT_PREFIX}.eigenvec"
    echo "  - Eigenvalues:  ${OUT_PREFIX}.eigenval"
else
    echo "PCA failed: Output files not found."
    exit 1
fi
