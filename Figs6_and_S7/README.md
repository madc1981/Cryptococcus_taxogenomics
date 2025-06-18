# rRNA Gene Analysis Pipeline

A unified Python-based pipeline for detecting and visualizing rRNA genes across fungal genomes, with statistical analyses for genomic distribution patterns.

## Features

- Annotates 5S, 5.8S, 18S, and 28S rRNA genes using [Barrnap](https://github.com/tseemann/barrnap)
- Detects canonical rDNA clusters (28S–5.8S–18S ± 5S)
- Generates genome-wide and zoomed-in plots of rRNA gene positions
- Tests for non-random chromosomal distribution of 5S genes
- Assesses 5S subtelomeric enrichment using Fisher’s exact test
  - Subtelomeric regions can be defined by fixed length (bp) or proportional to chromosome size
- Supports optional centromere annotation and custom plotting regions
- Outputs structured TSV tables and SVG figures per genome

## Requirements

- Python ≥ 3.8
- External tools:
  - [Barrnap](https://github.com/tseemann/barrnap)
  - Samtools
- Python packages:
  - `pandas`
  - `matplotlib`
  - `scipy`
  - `cairosvg`

### Recommended installation via Conda

```bash
conda create -n rRNA_pipeline_env -c conda-forge -c bioconda \
    python=3.10 samtools barrnap pandas matplotlib cairosvg scipy
conda activate rRNA_pipeline_env
```

## Usage

```bash
python 0_rRNA_pipeline.py --input_dir genomes/ --threads 16 \
    [--centromere_dir centromeres/] [--defined_order order.txt] \
    [--custom_regions_file zooms.tsv] [--output_dir results/] \
    [--subtelomere_window 50000] [--subtelomere_percent 0.025]
```

### Input Requirements

- Genome files in FASTA format, named as `*.scaffolds.fa` or `*.contigs.fa`
- Optional centromere files in BED format: `<strain>.cen.bed`, `<strain>.centromeres.bed`, or `<strain>.centro.bed`
- Optional zoom region file (`TSV`) with four columns: `Strain`, `Contig`, `Start`, `End`

### Notes

- Use `--subtelomere_window` to define subtelomeric regions as a fixed number of base pairs from each contig end (default: 50,000 bp)
- Alternatively, use `--subtelomere_percent` to define subtelomeric regions as a proportion of each contig length (e.g., `0.025` = 2.5% per end)
- Only one of these two options should be used per run
- Run with `--help` to view all available command-line options

## Outputs

- Genome-wide rRNA distribution plots (`SVG`)
- Zoomed-in cluster plots (per-genome and combined PDF)
- TSV tables summarizing:
  - rRNA gene coordinates
  - rDNA cluster counts
  - 5S chromosomal distribution and enrichment
  - 5S subtelomeric enrichment results

## Citation

If you use this pipeline, please cite:
> Coelho et al. 2025. *Genomic and phenotypic insights into the expanding phylogenetic landscape of the Cryptococcus genus* (in preparation).

## Contact

For questions or contributions, please contact: marco.dias.coelho [at] duke.edu.
