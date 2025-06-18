#!/usr/bin/env python3
# ==============================================================================
# Author: Marco A. Coelho @ Heitman Lab, Duke University
# Date: 2025-03-31
# Version: 1.2 (publication-ready)
# Description:
#     Correlation analysis between ANI and AAI values across Cryptococcus strains.
#     Merges ANI and AAI matrices with metadata, performs correlation tests,
#     highlights selected species comparisons, and produces annotated scatter plots.
#
# Outputs:
#     - <prefix>_correlation_plot.pdf
#     - <prefix>_divergence_plot.pdf
#     - <prefix>_correlation_stats.txt
#     - <prefix>_top_interspecies_identity.csv
#
# Requirements:
#     - pandas
#     - seaborn
#     - matplotlib
#     - scipy
#
# Usage:
#     python 0_plot_ani_vs_aai_correlation.py \
#         --ani formatted_ani.tsv \
#         --aai computed_aai.tsv \
#         --metadata metadata.tsv \
#         --specific_comparisons C_baci,C_deca \
#         --prefix Cryptococcus
# ==============================================================================

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr
import argparse
import sys


def parse_args():
    parser = argparse.ArgumentParser(description="Compute correlation between ANI and AAI.")
    parser.add_argument('--ani', required=True, help='Path to ANI matrix (wide format).')
    parser.add_argument('--aai', required=True, help='Path to AAI table.')
    parser.add_argument('--metadata', required=True, help='Path to strain metadata file.')
    parser.add_argument('--specific_comparisons', help='Comma-separated list of species to highlight.', default=None)
    parser.add_argument('--prefix', default="ani_aai", help='Prefix for all output files')
    return parser.parse_args()


def load_and_process_ani(path):
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index.name = "Label 1"
    df = df.reset_index().melt(id_vars=["Label 1"], var_name="Label 2", value_name="ANI")
    df = df[df["Label 1"] != df["Label 2"]]
    df["key"] = df.apply(lambda x: tuple(sorted([x["Label 1"], x["Label 2"]])), axis=1)
    df = df.drop_duplicates("key").drop(columns=["key"])
    return df


def load_and_merge_data(ani_df, aai_path, metadata_path):
    aai_df = pd.read_csv(aai_path, sep="\t")
    merged = pd.merge(ani_df, aai_df[["Label 1", "Label 2", "AAI"]], on=["Label 1", "Label 2"], how="inner")

    metadata = pd.read_csv(metadata_path, sep="\t")
    merged = merged.merge(metadata[["Strain", "Species"]], left_on="Label 1", right_on="Strain", how="left").drop(columns=["Strain"])
    merged = merged.rename(columns={"Species": "Species_1"})
    merged = merged.merge(metadata[["Strain", "Species"]], left_on="Label 2", right_on="Strain", how="left").drop(columns=["Strain"])
    merged = merged.rename(columns={"Species": "Species_2"})

    merged["Comparison Type"] = merged["Species_1"] == merged["Species_2"]
    merged["Comparison Type"] = merged["Comparison Type"].map({True: "intraspecies", False: "interspecies"})

    return merged


def compute_and_save_correlations(df, prefix):
    pearson_r, pearson_p = pearsonr(df["ANI"], df["AAI"])
    spearman_rho, spearman_p = spearmanr(df["ANI"], df["AAI"])

    with open(f"{prefix}_correlation_stats.txt", "w") as f:
        f.write(f"Spearman correlation (\u03c1) = {spearman_rho:.4f}, p = {spearman_p:.4e}\n")
        f.write(f"Pearson correlation (r) = {pearson_r:.4f}, p = {pearson_p:.4e}\n")

    return pearson_r, pearson_p, spearman_rho, spearman_p


def annotate_specific_comparisons(df, species_list):
    if species_list is None:
        df["Highlight"] = "Other"
    else:
        species_set = set(species_list)
        def tag(row):
            sp1 = row["Species_1"]
            sp2 = row["Species_2"]
            if sp1 in species_set and sp2 in species_set and sp1 != sp2:
                return "Specific"
            return "Other"
        df["Highlight"] = df.apply(tag, axis=1)
    return df


def plot_identity(df, pearson_r, pearson_p, spearman_rho, spearman_p, prefix):
    plt.figure(figsize=(8, 6))
    style_val = "Highlight" if "Specific" in df["Highlight"].values else None
    size_val = "Highlight" if "Specific" in df["Highlight"].values else None
    size_map = {"Specific": 100, "Other": 20} if "Specific" in df["Highlight"].values else None

    sns.scatterplot(
        data=df,
        x="ANI", y="AAI",
        hue="Comparison Type",
        style=style_val,
        size=size_val,
        sizes=size_map,
        alpha=0.6,
        edgecolor=None
    )
    plt.title(
        f"ANI vs AAI by species comparison\n"
        f"Spearman ρ = {spearman_rho:.3f} (p = {spearman_p:.1e})   "
        f"Pearson r = {pearson_r:.3f} (p = {pearson_p:.1e})"
    )
    plt.xlabel("Average Nucleotide Identity (%)")
    plt.ylabel("Average Amino Acid Identity (%)")
    plt.grid(True)
    plt.legend(title="Comparison / Highlight")
    plt.tight_layout()
    plt.rcParams['pdf.fonttype'] = 42
    plt.savefig(f"{prefix}_ani_vs_aai_correlation_plot.pdf")
    plt.close()


def plot_divergence(df, pearson_r, pearson_p, spearman_rho, spearman_p, prefix):
    df["Divergence_ANI"] = 100 - df["ANI"]
    df["Divergence_AAI"] = 100 - df["AAI"]

    plt.figure(figsize=(8, 6))
    style_val = "Highlight" if "Specific" in df["Highlight"].values else None
    size_val = "Highlight" if "Specific" in df["Highlight"].values else None
    size_map = {"Specific": 100, "Other": 20} if "Specific" in df["Highlight"].values else None

    sns.scatterplot(
        data=df,
        x="Divergence_ANI", y="Divergence_AAI",
        hue="Comparison Type",
        style=style_val,
        size=size_val,
        sizes=size_map,
        alpha=0.6,
        edgecolor=None
    )
    plt.title(
        f"Genomic divergence: 100 - ANI vs 100 - AAI\n"
        f"Spearman ρ = {-spearman_rho:.3f} (p = {spearman_p:.1e})   "
        f"Pearson r = {-pearson_r:.3f} (p = {pearson_p:.1e})"
    )
    plt.xlabel("Nucleotide divergence (100 - ANI %)")
    plt.ylabel("Amino acid divergence (100 - AAI %)")
    plt.grid(True)
    plt.legend(title="Comparison / highlight")
    plt.tight_layout()
    plt.rcParams['pdf.fonttype'] = 42
    plt.savefig(f"{prefix}_ani_vs_aai_divergence_plot.pdf")
    plt.close()


def save_top_interspecies(df, prefix):
    df["Mean Identity"] = (df["ANI"] + df["AAI"]) / 2
    top = df[df["Comparison Type"] == "interspecies"]
    top = top.sort_values(by="Mean Identity", ascending=False).head(20)
    top.to_csv(f"{prefix}_top_interspecies_identity.csv", index=False)


def main():
    args = parse_args()

    ani_long = load_and_process_ani(args.ani)
    merged = load_and_merge_data(ani_long, args.aai, args.metadata)

    species_list = args.specific_comparisons.split(",") if args.specific_comparisons else None
    merged = annotate_specific_comparisons(merged, species_list)

    pearson_r, pearson_p, spearman_rho, spearman_p = compute_and_save_correlations(merged, args.prefix)
    plot_identity(merged, pearson_r, pearson_p, spearman_rho, spearman_p, args.prefix)
    plot_divergence(merged, pearson_r, pearson_p, spearman_rho, spearman_p, args.prefix)
    save_top_interspecies(merged, args.prefix)


if __name__ == "__main__":
    main()