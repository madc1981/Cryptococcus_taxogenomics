# Biolog Phenotype Analysis Pipeline

Custom Python pipeline for analyzing Biolog substrate utilization data across Cryptococcus strains

## Overview
The script `0_biolog_analysis_pipeline.py` processes a curated Excel matrix of phenotypic growth responses, performs dimensionality reduction (PCA), hierarchical clustering, and machine learning (Random Forest) to identify phenotypic traits informative of clade-level differentiation. It also generates summary plots and flags redundant traits.

## Input Format
- Input must be an `.xlsx` file (e.g., `Cryptococcus_Biolog_data_parsed.xlsx`) where:
  - Rows represent phenotypic assays (e.g., substrates).
  - Columns represent Cryptococcus strains.
  - The first column ("Substrate") contains test names.
  - Other cells should be annotated as "positive", "negative", or left blank.
  - Strain names must include a clade prefix (e.g., `A_Cneo_H99`), indicating clade A strain.

## Output
The script generates the following key outputs:
- PCA plot (`*_pca.pdf`) with 95% confidence ellipses by clade
- Clustermap (`*_phenotype_clustermap.pdf`) and optional heatmap
- Top 15 informative phenotypes by Random Forest (`*_rf_feature_importance.csv`)
- Barplot of clade-specific phenotype frequencies (`*_top15_phenotypes_by_clade_barplot.pdf/png`)
- List of redundant/duplicate phenotypes (`*_duplicate_phenotypes.csv`)

## Conda Environment
To ensure reproducibility, a Conda environment is recommended and the yml file is provided.
Create it with:

```bash
conda env create -f environment.yml
conda activate biolog_env
```

## Usage
```bash
python 0_biolog_analysis_pipeline.py Cryptococcus_Biolog_data_parsed.xlsx
```

## Citation

If you use this pipeline, please consider citing:
> Coelho et al. 2025. *Genomic and phenotypic insights into the expanding phylogenetic landscape of the Cryptococcus genus* (in preparation).

---
For questions or contributions, please contact: marco.dias.coelho [at] duke.edu
