#!/bin/bash

# ==============================================================================
# Author: Marco A. Coelho @ Heitman Lab, Duke University
# Date: 2025-03-27
# Version: 1.4
# Description: Pipeline to compute average amino acid identity (AAI) using EzAAI.
#              Converts .proteins.fa files to EzAAI-compatible databases and computes
#              pairwise AAI values with user-defined parameters.
#              Dendrogram generation is excluded to preserve interpretability of distances.
# Requirements: Conda environment 'ezaai_env' with EzAAI v1.2.3 and supported tools (mmseqs, diamond, or blastp)
# Usage: ./0_run_ezaai_pipeline.sh [PROGRAM] [IDENTITY] [COVERAGE] [PREFIX]
# Example: ./0_run_ezaai_pipeline.sh diamond 0.4 0.5 Cryptococcus
# ==============================================================================

# Default parameters
DEFAULT_PROGRAM="mmseqs"
DEFAULT_IDENTITY="0.3"
DEFAULT_COVERAGE="0.5"
THREADS=28

# Input parameters from command line (with fallbacks)
PROGRAM="${1:-$DEFAULT_PROGRAM}"
IDENTITY="${2:-$DEFAULT_IDENTITY}"
COVERAGE="${3:-$DEFAULT_COVERAGE}"
DEFAULT_PREFIX="aai"

# Extend the command-line argument handling
PREFIX="${4:-$DEFAULT_PREFIX}"

# Output filenames (include params for clarity)
AAI_OUTPUT="${PREFIX}_${PROGRAM}_id${IDENTITY}_cov${COVERAGE}.tsv"


echo "Running EzAAI with:"
echo "  Program:    $PROGRAM"
echo "  Identity:   $IDENTITY"
echo "  Coverage:   $COVERAGE"
echo "  Threads:    $THREADS"

# Create necessary directories
mkdir -p db tmp

# Activate conda environment
eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
conda activate ezaai_env || { echo "Failed to activate ezaai_env"; exit 1; }

echo "Starting database conversion..."

# Convert all .proteins.fa files to .db
for file in *.proteins.fa; do
    base_name=$(basename "$file" .proteins.fa)
    ezaai convert -i "$file" -s prot -o "db/${base_name}.db" -l "$base_name" -tmp tmp
    echo "Converted: $file -> db/${base_name}.db"
done

echo "Database conversion completed!"

# Run AAI calculation
echo "Running AAI calculation with $PROGRAM..."
ezaai calculate -i db/ -j db/ -o "$AAI_OUTPUT" -p "$PROGRAM" -t "$THREADS" -id "$IDENTITY" -cov "$COVERAGE" -self 1

echo "AAI calculation completed! Results saved in $AAI_OUTPUT"

# Deactivate environment
conda deactivate
