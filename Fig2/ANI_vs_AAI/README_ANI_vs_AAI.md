# ANI vs AAI Correlation Analysis

This folder contains the script and input requirements for generating correlation and divergence plots comparing average nucleotide identity (ANI) and average amino acid identity (AAI) across Cryptococcus genomes.

## Overview

The pipeline:

1. Merges pairwise ANI and AAI data with species metadata.
2. Calculates Pearson and Spearman correlation coefficients.
3. Highlights selected species/species comparisons (optional).
4. Generates annotated scatter plots for both identity and divergence metrics.

## File List

```
ANI_vs_AAI/
├── 0_plot_ani_vs_aai_correlation.py   # Script to analyze and visualize ANI vs AAI relationships
├── README.md
```

## Requirements

- Python 3.7+
- pandas
- seaborn
- matplotlib
- scipy

Install dependencies:

```bash
pip install pandas seaborn matplotlib scipy
```

## Inputs

- `--ani`: Tab-delimited ANI matrix in wide format (as formatted by ANI pipeline).
- `--aai`: Tab-delimited AAI table in long format (as produced by EzAAI).
- `--metadata`: Metadata table with at least three columns: `Strain`, `Clade`, and `Species` (used for grouping and optional highlights).
- `--specific_comparisons` *(optional)*: Comma-separated species names to highlight.
- `--prefix` *(optional)*: Prefix for output files (default: `ani_aai`).

## Example

```bash
python 0_plot_ani_vs_aai_correlation.py \
    --ani Cryptococcus_ani_matrix.tsv \
    --aai Cryptococcus_aai_diamond_id0.4_cov0.5.tsv \
    --metadata Cryptococcus_metadata.tsv \
    --specific_comparisons Cryptococcus_bacillisporus,Cryptococcus_decagattii \
    --prefix Cryptococcus
```

## Outputs

- `Cryptococcus_correlation_plot.pdf`: ANI vs AAI annotated scatter plot.
- `Cryptococcus_divergence_plot.pdf`: Divergence scatter plot (100 – ANI vs 100 – AAI).
- `Cryptococcus_correlation_stats.txt`: Pearson and Spearman correlation coefficients.
- `Cryptococcus_top_interspecies_identity.csv`: Top 20 interspecies comparisons based on combined ANI and AAI.

## Citation

If you use this script or analysis pipeline in your research, please cite:

- Coelho et al. 2025. *Genomic and phenotypic insights into the expanding phylogenetic landscape of the Cryptococcus genus* (in preparation).

---

For questions or issues, contact: marco.dias.coelho [at] duke.edu

