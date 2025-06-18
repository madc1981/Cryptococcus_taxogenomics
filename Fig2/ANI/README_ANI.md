# ANI Analysis Pipeline

This subfolder contains the scripts and instructions used to compute average nucleotide identity (ANI) and generate clustering dendrograms to support Figure 2 of the manuscript.

## Overview

This pipeline is designed to:

1. Parse OrthoANIu output matrices and reformat them.
2. Generate distance-based dendrograms and rooted Newick trees from ANI values.
3. Compute minimum, maximum, and mean ANI values between species or clade groupings.

## File List

```
ANI/
├── 0_parse_orthoANI_matrix_to_table.py       # Format OrthoANIu matrix to TSV
├── 1_plot_dendrogram_and_tree.py             # Generate dendrogram and Newick tree
├── 2_compute_species_pairwise_ani.py         # Compute min/max/mean ANI by group
├── README.md
```

## Requirements

- Python 3.7+
- pandas
- matplotlib
- scipy
- numpy

Install dependencies using pip:

```bash
pip install pandas matplotlib scipy numpy
```

## Usage

### 1. Format the OrthoANIu matrix

```bash
python 0_parse_orthoANI_matrix_to_table.py path/to/raw_matrix.txt > formatted_matrix.tsv
```

Optionally exclude a strain:

```bash
python 0_parse_orthoANI_matrix_to_table.py path/to/raw_matrix.txt CBS12345
```

### 2. Plot dendrogram and generate Newick tree

```bash
python 1_plot_dendrogram_and_tree.py formatted_matrix.tsv
```

This script computes a distance matrix by subtracting %ANI from 100, performs average linkage hierarchical clustering (UPGMA-like), and exports both a dendrogram (PDF) and rooted Newick-format tree.

Outputs:

- `formatted_matrix_dendrogram_v1a.pdf`
- `formatted_matrix_tree.nwk`

### 3. Compute ANI stats between species or clades

```bash
python 2_compute_species_pairwise_ani.py formatted_matrix.tsv strain_to_clade.tsv
```

This produces a CSV file summarizing minimum, maximum, and mean ANI values between all clade or species pairs defined in the metadata.

## Input Files

- **OrthoANIu matrix file**: Raw output from [https://www.ezbiocloud.net/tools/orthoaniu](https://www.ezbiocloud.net/tools/orthoaniu)
- **strain\_to\_clade.tsv** *(optional)*: Two-column tab-separated file mapping strain names to clades or species labels.

Example:

```
Cryptococcus_neoformans_H99	CladeVNI
Cryptococcus_gattii_WM276	Clade5
Cryptococcus_depauperatus_CBS7855	Clade14
```

## Citation

If you use OrthoANIu in your research, please cite:

- Lee I, Kim YO, Park SC, Chun J. *OrthoANI: An improved algorithm and software for calculating average nucleotide identity.* Int J Syst Evol Microbiol. 2016 Jan;66(1):1100–1103. doi: 10.1099/ijsem.0.000760. PMID: 26585518.

If you use this ANI pipeline and accompanying scripts, please cite:

- Coelho et al. 2025. *Genomic and phenotypic insights into the expanding phylogenetic landscape of the Cryptococcus genus* (in preparation).

## License

The scripts are released under the MIT License.

---

For questions or issues, please contact: marco.dias.coelho [at] duke.edu