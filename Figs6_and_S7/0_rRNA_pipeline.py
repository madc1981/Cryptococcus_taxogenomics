#!/usr/bin/env python3

# =============================================================================================================
# Author: Marco A. Coelho @ Heitman lab, Duke University
# Date: 2025-05-27
# Version: 3.1
#
# Description:
#   Unified rRNA analysis pipeline that:
#     - Detects rRNA genes using Barrnap
#     - Generates scaled genome-wide chromosome plots showing rRNA genes (5S, 5.8S, 18S, 28S)
#     - Identifies and plots canonical rDNA clusters (28S–5.8S–18S) ± 5S
#     - Tests 5S chromosomal distribution enrichment using chi-squared goodness-of-fit
#     - Tests subtelomeric enrichment of 5S genes using Fisher’s exact test
#         - Subtelomeric regions can be defined as either a fixed window (in bp) or a proportion of each contig
#     - Summarizes rRNA gene content and outputs per-genome and combined TSVs
#     - Supports optional karyotype ordering and user-defined zoom-in regions
#
# Usage:
#   python 0_rRNA_pipeline.py --input_dir genomes/ --threads 16 \
#       [--centromere_dir centromeres/] [--defined_order order.txt] \
#       [--custom_regions_file zooms.tsv] [--output_dir results/] \
#       [--subtelomere_window 50000] [--subtelomere_percent 0.025]
#
#   - All input genomes should be in FASTA format and named as *.scaffolds.fa or *.contigs.fa.
#   - Corresponding Barrnap outputs, rRNA coordinates, and centromere files (optional) will be inferred
#     from filenames. Output files are written per-genome in TSV and SVG format.
#
#   - If provided, centromere coordinates should be in BED format with three tab-separated columns:
#        contig_name    start_position    end_position
#     Centromere files must be named as *.cen.bed, and placed in the directory specified by --centromere_dir.
#
#   - If provided, custom zoom regions must be specified in a tab-delimited file with four columns:
#        Strain    Contig    Start    End
#     Each row defines a specific genomic window to highlight and plot. This file is passed using --custom_regions_file.
#
#   - Subtelomeric regions can be defined either as a fixed size in base pairs (using --subtelomere_window)
#     or proportionally (e.g., 0.025, representing 2.5% of each contig end) using --subtelomere_percent.
#     Only one should be specified. If both are omitted, the default is 50,000 bp windows per end.
#
#   - For a full list of arguments and options, run:
#       python 0_rRNA_pipeline.py --help
#
# Dependencies:
#   - Python >= 3.8
#   - External tools:
#       - Barrnap (https://github.com/tseemann/barrnap)
#       - Samtools
#   - Python libraries:
#       - pandas
#       - matplotlib
#       - cairosvg
#       - scipy
#
# Installation (recommended):
#   conda create -n rRNA_pipeline_env -c conda-forge -c bioconda \
#       python=3.10 samtools barrnap pandas matplotlib cairosvg scipy
#   conda activate rRNA_pipeline_env
# =============================================================================================================



import os
import glob
import subprocess
import argparse
import csv
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import Ellipse
from matplotlib.ticker import FuncFormatter
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.image as mpimg
import cairosvg
from scipy.stats import chisquare
from scipy.stats import fisher_exact

# Global style for SVG/PDF exports
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'


def get_contig_lengths(fasta_file):
    """
    Parse FASTA .fai index to get contig lengths. Indexes file if missing.
    """
    fai_file = fasta_file + ".fai"
    if not os.path.exists(fai_file):
        print(f"Indexing {fasta_file} with samtools faidx...")
        subprocess.run(["samtools", "faidx", fasta_file], check=True)
    lengths = {}
    with open(fai_file, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            contig, length = parts[0], int(parts[1])
            lengths[contig] = length
    return lengths


def run_barrnap(genome_file, output_gff, output_fasta, threads):
    """
    Run Barrnap on input genome FASTA, save GFF and extracted rRNA sequences.
    """
    print(f"Running barrnap on {genome_file}...")
    with open(output_gff, "w") as gff_out:
        subprocess.run([
            "barrnap",
            "--kingdom", "euk",
            "--threads", str(threads),
            "--outseq", output_fasta
        ], stdin=open(genome_file, "r"), stdout=gff_out, check=True)


def parse_barrnap_gff(gff_file):
    """
    Parse Barrnap GFF3 output to extract rRNA gene features.
    """
    rRNAs = []
    with open(gff_file, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            feature = parts[2]
            attributes = parts[8]
            if feature == "rRNA":
                if "18S" in attributes:
                    gene = "18S"
                elif "28S" in attributes:
                    gene = "28S"
                elif "5.8S" in attributes or "5_8S" in attributes:
                    gene = "5.8S"
                elif "5S" in attributes:
                    gene = "5S"
                else:
                    continue
                rRNAs.append({
                    "contig": parts[0],
                    "start": int(parts[3]),
                    "end": int(parts[4]),
                    "strand": parts[6],
                    "gene": gene
                })
    return pd.DataFrame(rRNAs)


def count_clusters_strict_partition(df, max_cluster_span=5000):
    """
    Count rDNA clusters of 28S–5.8S–18S, with or without 5S, within a local genomic span.
    """
    core_genes = {"28S", "5.8S", "18S"}
    count_core_only = 0
    count_full = 0

    for contig, subdf in df.groupby("contig"):
        for strand, strand_df in subdf.groupby("strand"):
            genes_df = strand_df.sort_values("start").reset_index(drop=True)
            used = set()
            i = 0
            while i < len(genes_df):
                if i in used:
                    i += 1
                    continue
                gene_i = genes_df.loc[i, "gene"]
                pos_i = genes_df.loc[i, "start"]
                j = i + 1
                found = {gene_i}
                used_indices = [i]
                cluster_found = False

                while j < len(genes_df):
                    if genes_df.loc[j, "start"] - pos_i > max_cluster_span:
                        break
                    g = genes_df.loc[j, "gene"]
                    found.add(g)
                    used_indices.append(j)
                    if core_genes.issubset(found):
                        if "5S" in found:
                            count_full += 1
                        else:
                            count_core_only += 1
                        used.update(used_indices)
                        i = max(used_indices) + 1
                        cluster_found = True
                        break
                    j += 1

                if not cluster_found:
                    i += 1

    return count_core_only, count_full


def test_5s_enrichment(chr_counts, chr_lengths, strain, output_dir):
    """
    Perform chi-squared test of 5S enrichment vs chromosome length distribution.
    Saves observed vs. expected table with standardized residuals.
    """
    matched = set(chr_counts) & set(chr_lengths)
    if len(matched) < 2:
        return "NA", "NA", "Only one contig with 5S"

    obs = [chr_counts[c] for c in matched]
    total_5s = sum(obs)
    total_len = sum(chr_lengths[c] for c in matched)
    exp = [(chr_lengths[c] / total_len) * total_5s for c in matched]
    chi2, pval = chisquare(f_obs=obs, f_exp=exp)
    conclusion = "Enriched" if pval < 0.05 else "Not enriched"

    rows = []
    for c in matched:
        obs_c = chr_counts[c]
        exp_c = (chr_lengths[c] / total_len) * total_5s
        residual = obs_c - exp_c
        std_residual = (residual / (exp_c**0.5)) if exp_c > 0 else 0
        rows.append({
            "Contig": c,
            "Contig_Length": chr_lengths[c],
            "Observed_5S": obs_c,
            "Expected_5S": round(exp_c, 2),
            "Raw_Residual": round(residual, 2),
            "Standardized_Residual": round(std_residual, 2)
        })

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(output_dir, f"{strain}_5S_expected_vs_observed.tsv"), sep="\t", index=False)
    return round(chi2, 3), round(pval, 6), conclusion


def find_centromere_file(strain, centromere_dir):
    """
    Find a centromere BED file matching the strain name, if available.
    """
    patterns = [
        f"{centromere_dir}/{strain}.cen.bed",
        f"{centromere_dir}/{strain}.centromeres.bed",
        f"{centromere_dir}/{strain}.centro.bed"
    ]
    for pattern in patterns:
        files = glob.glob(pattern)
        if files:
            return files[0]
    return None


def parse_centromeres(bed_file):
    """
    Parse a BED file to extract centromere coordinates.
    """
    cens = {}
    if not bed_file or not os.path.exists(bed_file):
        return cens
    with open(bed_file, "r") as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            contig, start, end = parts[0], int(parts[1]), int(parts[2])
            cens[contig] = (start, end)
    return cens


def load_custom_regions(tsv_path):
    """
    Load user-defined zoom regions for plotting.
    """
    return pd.read_csv(tsv_path, sep="\t", header=None, names=["Strain", "Contig", "Start", "End"])


def plot_zoomed_cluster(rRNA_df, contig_lengths, strain, output_dir):
    """
    Create a zoomed-in plot around the main rDNA cluster (18S–5.8S–28S), 
    including any rRNA genes (5S, 5.8S, 18S, 28S) inside the detected region.
    """
    df_cluster = rRNA_df[rRNA_df['gene'].isin(['18S', '5.8S', '28S'])]
    if df_cluster.empty:
        return

    best_contig = df_cluster['contig'].value_counts().idxmax()
    sub_cluster = df_cluster[df_cluster['contig'] == best_contig]

    # Define zoomed region (5 kb padding)
    start = max(sub_cluster['start'].min() - 5000, 0)
    end = min(sub_cluster['end'].max() + 5000, contig_lengths[best_contig])

    # Now, select all rRNA genes (including 5S!) within this window
    df_zoom = rRNA_df[
        (rRNA_df['contig'] == best_contig) &
        (rRNA_df['start'] <= end) &
        (rRNA_df['end'] >= start)
    ]

    fig, ax = plt.subplots(figsize=(12, 2.5), constrained_layout=True)
    fig.suptitle(f"{strain} – rDNA cluster", fontsize=12)

    rRNA_colors = {
        "5S": "#e11f26",    # Red
        "5.8S": "#9b59b6",  # Purple
        "18S": "#1a75bb",   # Blue
        "28S": "#31b7af"    # Green
    }

    for _, row in df_zoom.iterrows():
        gene_start = row["start"]
        gene_end = row["end"]
        color = rRNA_colors.get(row["gene"], "gray")
        ax.add_patch(patches.Rectangle((gene_start, 0.35), gene_end - gene_start, 0.3,
                                       facecolor=color, edgecolor='none', zorder=10))

    ax.set_xlim(start, end)
    ax.set_ylim(0.2, 1.0)
    ax.set_yticks([])
    ax.set_xlabel(f"Genomic position on {best_contig} (bp)", fontsize=11, labelpad=15)
    ax.tick_params(axis='x', labelsize=12)

    for spine in ['top', 'right', 'left']:
        ax.spines[spine].set_visible(False)
    ax.spines['bottom'].set_linewidth(1.2)

    legend_patches = [patches.Patch(color=color, label=gene) for gene, color in rRNA_colors.items()]
    ax.legend(handles=legend_patches, loc="upper right", fontsize=9)

    out_path = os.path.join(output_dir, f"{strain}_rRNA_zoomed.svg")
    plt.savefig(out_path)
    print(f"Zoomed rDNA cluster plot saved to: {out_path}")


def plot_custom_zoomed_region(rRNA_df, strain, contig, start, end, output_dir):
    """
    Plot user-specified zoom-in region showing all rRNA genes.
    """
    df_zoom = rRNA_df[
        (rRNA_df["contig"] == contig) &
        (rRNA_df["start"] <= end) &
        (rRNA_df["end"] >= start)
    ]

    if df_zoom.empty:
        print(f"No rRNA genes found in {strain} {contig}:{start}-{end}")
        return

    rRNA_colors = {
        "5S": "#e11f26", "5.8S": "#9b59b6",
        "18S": "#1a75bb", "28S": "#31b7af"
    }

    fig, ax = plt.subplots(figsize=(12, 2.5), constrained_layout=True)
    fig.suptitle(f"{strain} – {contig}:{start}-{end}", fontsize=12)

    for _, row in df_zoom.iterrows():
        gene_start = row["start"]
        gene_end = row["end"]
        color = rRNA_colors.get(row["gene"], "gray")
        ax.add_patch(patches.Rectangle((gene_start, 0.35), gene_end - gene_start, 0.3,
                                       facecolor=color, edgecolor='none', zorder=10))

    ax.set_xlim(start, end)
    ax.set_ylim(0.2, 1.0)
    ax.set_yticks([])
    region_span = end - start
    x_unit = 1000 if region_span >= 1000 else 1
    x_label_unit = "kb" if x_unit == 1000 else "bp"
    ax.set_xlabel(f"Genomic position on {contig} ({x_label_unit})", fontsize=11, labelpad=15)
    ax.set_xticks(ax.get_xticks())
    ax.tick_params(axis='x', labelsize=12)
    ax.ticklabel_format(style='plain', axis='x')

    for spine in ['top', 'right', 'left']:
        ax.spines[spine].set_visible(False)
    ax.spines['bottom'].set_linewidth(1.2)

    legend_patches = [patches.Patch(color=color, label=gene) for gene, color in rRNA_colors.items()]
    ax.legend(handles=legend_patches, loc="upper right", fontsize=9)

    safe_contig = contig.replace("|", "_").replace("/", "_")
    out_file = os.path.join(output_dir, f"{strain}_{safe_contig}_{start}_{end}_custom_zoomed.svg")
    plt.savefig(out_file)
    print(f"🔍 Custom zoomed plot saved to: {out_file}")


def plot_rRNA_genes(fasta_file, gff_file, centromere_bed=None, output_dir="plots"):
    """
    Create a genome-wide plot with chromosomes, rRNA genes, and centromeres.
    """
    contig_lengths = get_contig_lengths(fasta_file)
    rRNA_df = parse_barrnap_gff(gff_file)
    centromeres = parse_centromeres(centromere_bed)

    contigs = list(contig_lengths.keys())
    max_length = max(contig_lengths.values())
    buffer = int(max_length * 0.05)
    max_plot = max_length + buffer

    unit = 1_000_000 if max_length > 1_000_000 else 1_000
    ylabel = "Length (Mb)" if unit == 1_000_000 else "Length (kb)"

    genome_name = os.path.basename(fasta_file)
    for suffix in [".scaffolds.fa", ".scaffolds.fasta", ".contigs.fa", ".contigs.fasta", ".fa", ".fasta"]:
        if genome_name.endswith(suffix):
            genome_name = genome_name.replace(suffix, "")
            break

    fig, ax = plt.subplots(figsize=(10, 10), constrained_layout=True)
    fig.suptitle(f"{genome_name}", fontsize=14)

    bar_width = 0.4
    spacing_within = 0.7
    fixed_ellipse_height = 80000

    rRNA_colors = {
        "5S": "#e11f26", "5.8S": "#9b59b6",
        "18S": "#1a75bb", "28S": "#31b7af"
    }

    for i, contig in enumerate(contigs):
        x = i * spacing_within
        height = contig_lengths[contig]
        ax.add_patch(patches.Rectangle((x - bar_width/2, 0), bar_width, height,
                                       facecolor="white", edgecolor='black'))

        ax.text(x, -buffer * 0.2, str(i + 1), ha='center', va='bottom', fontsize=14)

        if contig in centromeres:
            start, end = centromeres[contig]
            ax.add_patch(patches.Rectangle((x - bar_width/2, start), bar_width, end - start,
                                           facecolor='black', alpha=0.85, edgecolor='none', zorder=10))
            cent_mid = start + (end - start) / 2
            ax.add_patch(Ellipse((x, cent_mid),
                                 width=bar_width * 0.6,
                                 height=fixed_ellipse_height,
                                 facecolor='white', edgecolor='black',
                                 lw=1.5, zorder=11))

    for _, row in rRNA_df.iterrows():
        if row["contig"] not in contigs or row["gene"] not in rRNA_colors:
            continue
        i = contigs.index(row["contig"])
        x = i * spacing_within
        midpoint = (row["start"] + row["end"]) // 2
        color = rRNA_colors[row["gene"]]
        ax.plot([x - bar_width/4, x + bar_width/4], [midpoint, midpoint], color=color, lw=2, zorder=12)

    ax.invert_yaxis()
    ax.set_xlim(-1, len(contigs) * spacing_within + 1)
    ax.set_ylim(max_plot, -buffer * 1.0)
    ax.set_ylabel(ylabel, fontsize=18)
    ax.tick_params(axis='y', labelsize=14)
    ax.set_xticks([])
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x/unit:.1f}"))
    ax.yaxis.grid(True, linestyle='--', alpha=0.5)

    for spine in ["top", "right", "bottom"]:
        ax.spines[spine].set_visible(False)

    legend_patches = [patches.Patch(color=color, label=gene) for gene, color in rRNA_colors.items()]
    ax.legend(handles=legend_patches, loc="lower right", fontsize=14)

    output_file = os.path.join(output_dir, f"{genome_name}_rRNA_plot.svg")
    plt.savefig(output_file)
    print(f"Full genome rRNA plot saved to: {output_file}")



def plot_combined_genome_panels(genome_infos, output_path):
    """
    Create a combined panel showing multiple genomes horizontally with shared Y-axis.
    """
    fig_width = 2 + sum(len(info["contig_lengths"]) for info in genome_infos) * 0.8
    fig, ax = plt.subplots(figsize=(fig_width, 10), constrained_layout=True)

    spacing_within = 0.5
    spacing_between = 1.5
    bar_width = spacing_within * 0.5
    fixed_ellipse_height = 80000

    rRNA_colors = {
        "5S": "#e11f26",
        "5.8S": "#9b59b6",
        "18S": "#1a75bb",
        "28S": "#31b7af"
    }

    # Determine global maximum chromosome length
    global_max = max(
        max(info["contig_lengths"].values()) for info in genome_infos
    )
    buffer = int(global_max * 0.05)
    max_plot = global_max + buffer
    unit = 1_000_000 if global_max > 1_000_000 else 1_000
    ylabel = "Length (Mb)" if unit == 1_000_000 else "Length (kb)"

    x_offset = 0
    for idx, genome in enumerate(genome_infos):
        contigs = list(genome["contig_lengths"].keys())
        strain = genome["strain"]
        contig_lengths = genome["contig_lengths"]
        rRNA_df = genome["rRNA_df"]
        centromeres = genome["centromeres"]

        for i, contig in enumerate(contigs):
            x = x_offset + i * spacing_within
            height = contig_lengths[contig]
            ax.add_patch(patches.Rectangle((x - bar_width/2, 0), bar_width, height,
                                           facecolor="white", edgecolor="black"))

            # Contig label (numbered)
            ax.text(x, -buffer * 0.2, str(i + 1), ha='center', va='bottom', fontsize=14)

            if contig in centromeres:
                start, end = centromeres[contig]
                ax.add_patch(patches.Rectangle((x - bar_width/2, start), bar_width, end - start,
                                               facecolor='black', alpha=0.85, edgecolor='none', zorder=10))
                cent_mid = start + (end - start) / 2
                ax.add_patch(Ellipse((x, cent_mid),
                                     width=bar_width * 0.6,
                                     height=fixed_ellipse_height,
                                     facecolor='white', edgecolor='black',
                                     lw=1.5, zorder=11))

        for _, row in rRNA_df.iterrows():
            if row["contig"] not in contigs or row["gene"] not in rRNA_colors:
                continue
            i = contigs.index(row["contig"])
            x = x_offset + i * spacing_within
            midpoint = (row["start"] + row["end"]) // 2
            color = rRNA_colors[row["gene"]]
            ax.plot([x - bar_width/4, x + bar_width/4], [midpoint, midpoint],
                    color=color, lw=2, zorder=12)

        # Genome label centered above block
        block_center = x_offset + (len(contigs) - 1) * spacing_within / 2
        ax.text(block_center, -buffer * 1.2, strain, ha='center', va='bottom', fontsize=18, fontweight='normal')

        x_offset += len(contigs) * spacing_within + spacing_between

    # Styling
    ax.invert_yaxis()
    ax.set_xlim(-1, x_offset)
    ax.set_ylim(max_plot, -buffer * 1.2)
    ax.set_ylabel(ylabel, fontsize=18)
    ax.tick_params(axis='y', labelsize=14)
    ax.set_xticks([])
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x/unit:.1f}"))
    ax.yaxis.grid(True, linestyle='--', alpha=0.5)

    for spine in ["top", "right", "bottom"]:
        ax.spines[spine].set_visible(False)

    legend_patches = [patches.Patch(color=color, label=gene) for gene, color in rRNA_colors.items()]
    ax.legend(handles=legend_patches, loc="lower right", fontsize=10)


    # Create mapping and export
    mapping_path = os.path.join(os.path.dirname(output_path), "contig_number_mapping.csv")
    with open(mapping_path, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Genome", "Original_Contig", "Assigned_Number"])
        for genome in genome_infos:
            strain = genome["strain"]
            contigs = list(genome["contig_lengths"].keys())
            for idx, contig in enumerate(contigs):
                writer.writerow([strain, contig, idx + 1])

    print(f"Contig mapping table saved to: {mapping_path}")

    plt.savefig(output_path)
    print(f"Combined full-genome plot saved to: {output_path}")


def perform_5s_enrichment_and_summary(rRNA_df, contig_lengths, strain, output_dir):
    """
    Count 5S per contig, perform enrichment test, and summarize all rRNA counts.
    """
    chr_counts = rRNA_df[rRNA_df["gene"] == "5S"]["contig"].value_counts().to_dict()
    chr_lengths = contig_lengths

    chi2, pval, conclusion = test_5s_enrichment(chr_counts, chr_lengths, strain, output_dir)
    df_chrom_dist = pd.DataFrame([
        {"Strain": strain, "Contig": contig, "5S_Count": count}
        for contig, count in chr_counts.items()
    ])
    df_chrom_dist.to_csv(os.path.join(output_dir, f"{strain}_5S_chromosomal_distribution.tsv"),
                         sep="\t", index=False)

    df_stats = pd.DataFrame([{
        "Strain": strain, "Chi2": chi2, "p_value": pval, "Conclusion": conclusion
    }])
    df_stats.to_csv(os.path.join(output_dir, f"{strain}_5S_enrichment_statistics.tsv"),
                    sep="\t", index=False)

    # rRNA summary counts
    gene_counts = rRNA_df["gene"].value_counts().to_dict()
    core, full = count_clusters_strict_partition(rRNA_df)
    gene_counts["rDNA_Cluster_28S-5.8S-18S"] = core
    gene_counts["rDNA_Cluster_28S-5.8S-18S-5S"] = full

    df_summary = pd.DataFrame([
        {"Strain": strain, "Gene": gene, "Count": gene_counts[gene]}
        for gene in sorted(gene_counts)
    ])
    df_summary.to_csv(os.path.join(output_dir, f"{strain}_rRNA_summary_counts.tsv"),
                      sep="\t", index=False)

def test_5s_subtelomeric_enrichment(rRNA_df, contig_lengths, strain, output_dir, window_size=None, window_percent=None):
    """
    Test if 5S rRNA genes are enriched in subtelomeric regions using Fisher’s exact test.
    Either a fixed window size (in bp) or a window_percent (proportion of chromosome length per end) must be provided.
    """
    if (window_size is None and window_percent is None) or (window_size and window_percent):
        raise ValueError("Specify either --subtelomere_window (bp) or --subtelomere_percent (float), not both.")

    total_genome_length = sum(contig_lengths.values())
    total_subtelomere_length = 0
    total_5s = len(rRNA_df[rRNA_df["gene"] == "5S"])
    subtelomeric_5s = 0

    for contig, length in contig_lengths.items():
        if window_percent is not None:
            flank_size = int(length * window_percent)
        else:
            flank_size = window_size

        subtelomere_regions = [
            (0, min(flank_size, length)),
            (max(0, length - flank_size), length)
        ]
        total_subtelomere_length += sum(end - start for start, end in subtelomere_regions)

        df_contig = rRNA_df[(rRNA_df["contig"] == contig) & (rRNA_df["gene"] == "5S")]
        for _, row in df_contig.iterrows():
            mid = (row["start"] + row["end"]) // 2
            for start, end in subtelomere_regions:
                if start <= mid <= end:
                    subtelomeric_5s += 1
                    break

    non_subtelomeric_5s = total_5s - subtelomeric_5s
    subtelomere_bp = total_subtelomere_length
    non_subtelomere_bp = total_genome_length - total_subtelomere_length

    contingency = [
        [subtelomeric_5s, non_subtelomeric_5s],
        [subtelomere_bp, non_subtelomere_bp]
    ]

    _, p_value = fisher_exact(contingency, alternative="greater")
    enriched = "Enriched" if p_value < 0.05 else "Not enriched"

    output_file = os.path.join(output_dir, f"{strain}_5S_subtelomeric_enrichment.tsv")
    df_result = pd.DataFrame([{
        "Strain": strain,
        "5S_Subtelomeric": subtelomeric_5s,
        "5S_Total": total_5s,
        "Subtelomeric_bp": subtelomere_bp,
        "Genome_bp": total_genome_length,
        "p_value": f"{p_value:.2e}",
        "Conclusion": enriched
    }])
    df_result.to_csv(output_file, sep="\t", index=False)
    print(f"Subtelomeric 5S enrichment results saved to: {output_file}")


def main():
    parser = argparse.ArgumentParser(description="Detect and plot rRNA genes from genome FASTA files.")
    parser.add_argument("--input_dir", type=str, required=True, help="Directory containing genome FASTA files")
    parser.add_argument("--threads", type=int, required=True, help="Number of threads for Barrnap")
    parser.add_argument("--centromere_dir", type=str, default=None, help="Directory containing centromere BED files (optional)")
    parser.add_argument("--defined_order", type=str, default=None, help="Optional file listing genome names (one per line) to define order in combined plot")
    parser.add_argument("--custom_regions_file", type=str, default=None, help="Optional TSV file with custom zoom regions: Strain,Contig,Start,End")
    parser.add_argument("--output_dir", type=str, default=None, help="Optional output directory (default = input_dir)")
    parser.add_argument("--subtelomere_window", type=int, default=50000,
                    help="Size of subtelomeric region from each chromosome end (default: 50,000 bp)")
    parser.add_argument("--subtelomere_percent", type=float, default=None,
    help="Proportion of each chromosome length to define as subtelomeric per end (e.g., 0.05 = 5%). Overrides --subtelomere_window if set.")

    args = parser.parse_args()


    input_dir = args.input_dir
    output_dir = args.output_dir if args.output_dir else input_dir
    centromere_dir = args.centromere_dir
    threads = args.threads
    defined_order_file = args.defined_order

    fasta_extensions = (".scaffolds.fa", ".scaffolds.fasta", ".contigs.fa", ".contigs.fasta")
    fasta_files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.endswith(fasta_extensions)]

    if not fasta_files:
        print(f"No genome FASTA files found in {input_dir}")
        return

    os.makedirs(output_dir, exist_ok=True)

    genome_infos = []
    zoomed_plots = []

    for fasta in fasta_files:
        base = os.path.basename(fasta)
        strain = base.replace(".scaffolds.fa", "").replace(".scaffolds.fasta", "") \
                     .replace(".contigs.fa", "").replace(".contigs.fasta", "")
        gff = os.path.join(input_dir, f"{strain}.rRNA.gff")
        rrna_fasta = os.path.join(input_dir, f"{strain}.rRNA.fa")
        cen_bed = find_centromere_file(strain, centromere_dir) if centromere_dir else None

        if not os.path.exists(gff):
            try:
                run_barrnap(fasta, gff, rrna_fasta, threads)
            except subprocess.CalledProcessError:
                print(f"Barrnap failed for {fasta}. Skipping.")
                continue

        contig_lengths = get_contig_lengths(fasta)
        rRNA_df = parse_barrnap_gff(gff)

        # Zoomed custom regions
        if args.custom_regions_file:
            custom_regions_df = load_custom_regions(args.custom_regions_file)
            for _, region in custom_regions_df[custom_regions_df["Strain"] == strain].iterrows():
                plot_custom_zoomed_region(
                    rRNA_df=rRNA_df,
                    strain=strain,
                    contig=region["Contig"],
                    start=int(region["Start"]),
                    end=int(region["End"]),
                    output_dir=output_dir
                )

        centromeres = parse_centromeres(cen_bed) if cen_bed else {}

        # Save gene coordinates
        genes_output = os.path.join(output_dir, f"{strain}_rRNA_genes.tsv")
        rRNA_df.to_csv(genes_output, sep="\t", index=False)
        print(f"rRNA gene coordinates saved to: {genes_output}")

        # Per-genome plots
        plot_rRNA_genes(fasta, gff, cen_bed, output_dir)
        plot_zoomed_cluster(rRNA_df, contig_lengths, strain, output_dir)

        zoom_path = os.path.join(output_dir, f"{strain}_rRNA_zoomed.svg")
        if os.path.exists(zoom_path):
            zoomed_plots.append(zoom_path)

        # Store for combined panel
        genome_infos.append({
            "strain": strain,
            "contig_lengths": contig_lengths,
            "rRNA_df": rRNA_df,
            "centromeres": centromeres
        })

        # Run 5S enrichment + summary stats
        perform_5s_enrichment_and_summary(rRNA_df, contig_lengths, strain, output_dir)

        test_5s_subtelomeric_enrichment(
            rRNA_df=rRNA_df,
            contig_lengths=contig_lengths,
            strain=strain,
            output_dir=output_dir,
            window_size=args.subtelomere_window if args.subtelomere_percent is None else None,
            window_percent=args.subtelomere_percent
        )

    # Reorder genomes for combined plot if specified
    if defined_order_file:
        with open(defined_order_file) as f:
            desired_order = [line.strip() for line in f if line.strip()]
        genome_infos.sort(key=lambda g: desired_order.index(g["strain"]) if g["strain"] in desired_order else float('inf'))

    # Combined genome-wide panel
    if len(genome_infos) > 1:
        combined_output = os.path.join(output_dir, "combined_rRNA_all_genomes.svg")
        plot_combined_genome_panels(genome_infos, combined_output)

    # Combined zoomed-in panel (optional)
    if len(zoomed_plots) > 1:
        pdf_path = os.path.join(output_dir, "combined_rRNA_zoomed_plots.pdf")
        with PdfPages(pdf_path) as pdf:
            for svg_file in zoomed_plots:
                png_file = svg_file.replace(".svg", ".png")
                cairosvg.svg2png(url=svg_file, write_to=png_file)
                img = mpimg.imread(png_file)
                fig, ax = plt.subplots(figsize=(12, 2.5))
                ax.imshow(img)
                ax.axis('off')
                pdf.savefig(fig)
                plt.close()
                os.remove(png_file)
        print(f"Combined zoomed-in plots saved to: {pdf_path}")

    # Combine all output tables across strains 
    categories = [
        ("_5S_chromosomal_distribution.tsv", "all_strains_5S_chromosomal_distribution.tsv"),
        ("_5S_enrichment_statistics.tsv", "all_strains_5S_enrichment_statistics.tsv"),
        ("_5S_expected_vs_observed.tsv", "all_strains_5S_expected_vs_observed.tsv"),
        ("_rRNA_summary_counts.tsv", "all_strains_rRNA_summary_counts.tsv"),
        ("_5S_subtelomeric_enrichment.tsv", "all_strains_5S_subtelomeric_enrichment.tsv"),
    ]

    for suffix, merged_name in categories:
        all_tables = glob.glob(os.path.join(output_dir, f"*{suffix}"))
        dfs = []
        for file in all_tables:
            try:
                df = pd.read_csv(file, sep="\t")
                if suffix == "_5S_expected_vs_observed.tsv":
                    strain = os.path.basename(file).replace(suffix, "")
                    df.insert(0, "Strain", strain)
                dfs.append(df)
            except Exception as e:
                print(f"Could not read {file}: {e}")
        if dfs:
            merged_df = pd.concat(dfs, ignore_index=True)
            merged_df.to_csv(os.path.join(output_dir, merged_name), sep="\t", index=False)
            print(f"Merged table saved: {merged_name}")

if __name__ == "__main__":
    main()
