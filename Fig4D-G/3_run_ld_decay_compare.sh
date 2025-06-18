#!/bin/bash

#===============================================================================
# Author: Marco A. Coelho @ Heitman Lab, Duke University
# Date: 2025-05-28
# Version: 1.0
#
# Description:
#   Runs PopLDdecay on a multi-sample VCF file from Snippy with and without
#   minor allele frequency (MAF) filtering to generate LD decay statistics.
#
# Requirements:
#   - Conda environment: popgen_env
#   - Tool: PopLDdecay (https://github.com/BGI-shenzhen/PopLDdecay)
#
# Usage:
#   conda activate popgen_env
#   bash 3_run_ld_decay_compare.sh
#
# Input:
#   core.vcf                         # Multi-sample VCF
#
# Output:
#   core_lddecay_maf0.167_dist50000.stat.gz
#   core_lddecay_nomaf_dist50000.stat.gz
#===============================================================================

set -euo pipefail

# --- PARAMETERS ---
VCF_FILE="core.vcf"
OUT_PREFIX="core_lddecay"
MAXDIST=50000
MAF_FILTER=0.167

# --- CHECK TOOL ---
if ! command -v PopLDdecay &> /dev/null; then
    echo "Error: PopLDdecay not found in PATH. Please activate popgen_env."
    exit 1
fi

# --- RUN: With MAF filter ---
echo "LD decay with MAF ≥ $MAF_FILTER (MaxDist = $MAXDIST)..."
PopLDdecay -InVCF "$VCF_FILE" \
           -OutStat "${OUT_PREFIX}_maf${MAF_FILTER}_dist${MAXDIST}" \
           -MaxDist "$MAXDIST" \
           -MAF "$MAF_FILTER"

# --- RUN: No MAF filter ---
echo "LD decay without MAF filter (MaxDist = $MAXDIST)..."
PopLDdecay -InVCF "$VCF_FILE" \
           -OutStat "${OUT_PREFIX}_nomaf_dist${MAXDIST}" \
           -MaxDist "$MAXDIST"

# --- CHECK OUTPUT ---
if ls ${OUT_PREFIX}_maf${MAF_FILTER}_dist${MAXDIST}.stat* &>/dev/null && \
   ls ${OUT_PREFIX}_nomaf_dist${MAXDIST}.stat* &>/dev/null; then
    echo "LD decay completed:"
    echo "   → ${OUT_PREFIX}_maf${MAF_FILTER}_dist${MAXDIST}.stat.gz"
    echo "   → ${OUT_PREFIX}_nomaf_dist${MAXDIST}.stat.gz"
else
    echo "LD decay run failed: expected output files not found."
    exit 1
fi
