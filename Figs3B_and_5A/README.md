# Synteny analysis pipeline - 0_synteny_karyoplotter.py

## Overview

This script performs nucleotide-based synteny analysis between genome assemblies and produces high-quality chromosome-scale visualizations. It uses minimap2 to identify synteny blocks and generates scalable SVG plots that highlight conserved regions and optionally annotated centromeres. This script was used to generate raw Figures 3B and 5A of the associated manuscript (see citation below)

Visualizations include:
- To-scale chromosome representations
- Color-coded synteny blocks
- Centromere locations shown as black bars and white ellipses
- Optional multi-genome comparisons against a fixed reference

---

## Features

- Supports pairwise and reference-based synteny comparisons
- Annotates and visualizes centromere positions (if available)
- Outputs scalable SVG plots suitable for publication
- Handles multiple input formats and scales figure layout dynamically

---

## Requirements

- Python 3
- [minimap2](https://github.com/lh3/minimap2) (must be in your `$PATH`)
- Python libraries:
  - Biopython
  - pandas
  - matplotlib
  - numpy

Install dependencies using:
```bash
pip install biopython pandas matplotlib numpy
```

---

## Input Files

### Required

- Genome assemblies in FASTA format:
  - Supported extensions: `.fa`, `.fasta`
  - Each contig or chromosome must have a unique header.

Example:
```fasta
>chr_1
ATCGATCG...
>chr_2
GCTAGCTA...
```

### Optional (for centromere annotation)

- BED files for centromere coordinates (BED4 format):
  - Must be named with the same prefix as the genome file:
    ```
    Cryptococcus_neoformans_125.91.scaffolds.fa  →  Cryptococcus_neoformans_125.91.cen.bed
    ```

Example `cen.bed`:
```
chr_1   960921   1030510  CEN1
chr_2   904508   935491   CEN2
...
```

---

## Usage

### Basic pairwise comparison (all-vs-all):
```bash
python3 0_synteny_karyoplotter.py --input_dir /path/to/genomes
```

### Fixed-reference comparison (all vs one):
```bash
python3 0_synteny_karyoplotter.py --input_dir /path/to/genomes --ref Cryptococcus_neoformans_125.91
```

> The script uses `minimap2 -x asm20`, which is designed for genome pairs with up to ~10% sequence divergence. In practice, it performs well with higher divergence levels in fungal genomes.

---

## Output

- **Synteny block tables** (`*.tsv`)
- **SVG plots** for each comparison or multi-genome layout
- Output folders:
  - `results/` → minimap2 PAF files and parsed synteny tables
  - `plots/`   → final SVG visualizations

---

## Citation

If you use this pipeline, please cite:
> Coelho et al. 2025. *Genomic and phenotypic insights into the expanding phylogenetic landscape of the Cryptococcus genus* (in preparation).

---

For questions or issues, please contact: marco.dias.coelho [at] duke.edu

