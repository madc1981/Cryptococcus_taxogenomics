# Fig 4 - Pop Gen analysis of Cryptococcus porticicola

This folder contains scripts used to generate Figure 4 panels E, F and G of our manuscript.

All scripts were run on multi-sample VCF files generated by the **Snippy pipeline**, and results include:

- Principal Component Analysis (PCA)
- LD decay analysis with and without MAF filtering
- Four-gamete test to infer recombination signatures

---

## Requirements

All scripts run within a single Conda environment:

```bash
conda env create -f environment.yml
conda activate popgen_env
```

Tested with Python 3.11.

---

## Workflow Overview

| Step | Script                      | Purpose                                                    |
| ---- | --------------------------- | ---------------------------------------------------------- |
| 0    | `0_convert_vcf_to_plink.sh` | Converts `core.vcf` to PLINK binary format                 |
| 1    | `1_run_pca.sh`              | Runs PCA on `.bed/.bim/.fam` files with PLINK              |
| 2    | `2_plot_pca.py`             | Generates PC1 vs PC2 scatterplot                           |
| 3    | `3_run_ld_decay_compare.sh` | Computes LD decay with/without MAF filter using PopLDdecay |
| 4    | `4_plot_ld_decay.py`        | Plots LD decay curves, LOWESS smoothing, summary stats     |
| 5    | `5_four_gamete_test.py`     | Performs four-gamete test and barplot summary              |

All input/output filenames are hardcoded to match the `core.*` naming convention from Snippy pipeline.

---

## Inputs and Outputs

**Expected Input**:

- `core.vcf` — multi-sample VCF with high-quality biallelic SNPs from Snippy

**Main Outputs**:

- PCA plots: `core_pca_PC1_vs_PC2.pdf/.png`
- LD decay plots: `ld_decay_*.pdf/.png`, `ld_decay_summary.txt`
- Four-gamete plot: `four_gamete_test_plot.pdf`
- Intermediate stats: `.eigenvec`, `.stat.gz`, `.tsv`

---

## Citation

If you use or adapt this workflow, please cite the associated manuscript and tools:

- Snippy: [https://github.com/tseemann/snippy](https://github.com/tseemann/snippy)
- PLINK v1.9: [https://www.cog-genomics.org/plink/](https://www.cog-genomics.org/plink/)
- PopLDdecay: [https://github.com/BGI-shenzhen/PopLDdecay](https://github.com/BGI-shenzhen/PopLDdecay)
- VCFtools, pandas, seaborn, matplotlib, statsmodels

and, consider citing:
> Coelho et al. 2025. *Genomic and phenotypic insights into the expanding phylogenetic landscape of the Cryptococcus genus* (in preparation).

---

For questions or issues, please contact: marco.dias.coelho [at] duke.edu

