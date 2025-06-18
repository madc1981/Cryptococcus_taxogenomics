# Figure 2: ANI, AAI, dDDH, and Correlation Analyses

This folder contains all computational scripts, outputs, and readme files used for generating the ANI, AAI, dDDH, and ANI-vs-AAI correlation analyses featured in **Figure 2** of the manuscript.

Each subfolder contains:
- Scripts to run the analysis
- Supporting input and output files (processed data, plots, summary tables)
- A dedicated `README.md` with usage details and citations

---

## Folder Structure

```
Fig2/
├── AAI/
│   ├── 0_run_ezaii_pipeline.sh               # Main pipeline to compute AAI matrix
│   ├── 1_plot_aai_dendrogram.py              # Generates dendrogram and Newick tree from AAI values
│   ├── 2_compute_aai_stats.py                # Computes min/max/mean AAI values by clade/species
│   └── README_AAI.md
│
├── ANI/
│   ├── 0_parse_orthoANI_matrix_to_table.py   # Formats OrthoANIu matrix for analysis
│   ├── 1_plot_dendrogram_and_tree.py         # Creates ANI dendrogram and tree
│   ├── 2_compute_species_pairwise_ani.py     # Computes pairwise ANI stats
│   └── README_ANI.md
│
├── ANI_vs_AAI/
│   ├── 0_plot_ani_vs_aai_correlation.py      # Correlation plots and interspecies comparisons
│   └── README_ANI_vs_AAI.md
│
├── dDDH/
│   ├── 0_plot_dDDH_outputs_from_GGDC.py      # Summary and plots for GGDC Formula 2 dDDH results
│   └── README_dDDH.md
```

---

## Notes
- All scripts are standalone and parameterized.
- Each subfolder README provides versioning, input/output file formats, example usage, and appropriate citations.
- Scripts were executed independently; output files (e.g. plots, tables) are not listed here but are present within each folder.

---

## Dependencies

Each subfolder contains installation instructions and version requirements for the relevant tools.

## License

All scripts and documentation are released under the MIT License.

---

For questions or reproducibility inquiries, contact: marco.dias.coelho [at] duke.edu
