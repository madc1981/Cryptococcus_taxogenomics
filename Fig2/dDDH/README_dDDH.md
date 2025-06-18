# dDDH Analysis - GGDC results parser and visualization

This folder contains scripts and outputs for the digital DNA-DNA hybridization (dDDH) analysis of *Cryptococcus* genomes using Formula 2 results from the Genome-to-Genome Distance Calculator (GGDC) platform.

## Overview

Digital DNA-DNA hybridization (dDDH) values were computed using the GGDC web server ([https://ggdc.dsmz.de/ggdc.php](https://ggdc.dsmz.de/ggdc.php)), selecting specific intra- and inter-species comparisons. This analysis focuses on pairwise comparisons using GGDC's recommended Formula 2.

The script `0_plot_dDDH_outputs_from_GGDC.py` parses the GGDC output CSV, extracts strain and species labels, computes summary statistics, and generates selected visualizations.

## Input format

- The input **must be a single, compiled *.csv file** that merges all GGDC result tables into one.
- The file should preserve the original GGDC format, with column headers starting from the second row.
- Results must be obtained from [https://ggdc.dsmz.de/ggdc.php](https://ggdc.dsmz.de/ggdc.php).

## Output

The script produces the following outputs (based on the input file prefix):

### Tables:

- `*_formula2_summary_strain_labels.csv`: Clean summary table with DDH values, probabilities, and classification
- `*_inter_species_pair_legend.csv`: Annotated table with inter-species pair IDs for reference

### Plots:

- `*_ddh_formula2_barplot_strain_labels.pdf`: Bar plot of pairwise DDH values with intra- vs inter-species coloring
- `*_ddh_vs_probability_scatter_strain_labels.pdf`: Scatter plot of DDH (%) vs. probability DDH ≥ 70% (colored by comparison type)
- `*_scatter_intra_by_species.pdf`: Scatter plot of intra-species comparisons (colored by species)
- `*_scatter_inter_numbered.pdf`: Numbered inter-species scatter plot with pair labels

## Usage

Run the script on a single GGDC results CSV:

```bash
python3 0_plot_dDDH_outputs_from_GGDC.py Cryptococcus_ggdc_results_compiled.csv
```

This will generate the summary tables and the four selected plots described above. Only Formula 2 DDH values are considered.

## Requirements

- Python >= 3.7
- pandas
- seaborn
- matplotlib

## Citation

If you use GGDC in your research, please cite:

- Meier-Kolthoff JP, Sardà Carbasse J, Peinado-Olarte RL, Göker M. **TYGS and LPSN: a database tandem for fast and reliable genome-based classification and nomenclature of prokaryotes.** *Nucleic Acid Res* 2022 50:D801–D807.


If you use this analysis pipeline and associated scripts consider citing:

- Coelho et al. 2025. *Genomic and phenotypic insights into the expanding phylogenetic landscape of the Cryptococcus genus* (in preparation).

---

For questions, contact: marco.dias.coelho [at] duke.edu


