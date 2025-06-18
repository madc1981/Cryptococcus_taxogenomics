
#!/usr/bin/env python3

# ==============================================================================
# Script:     0_synteny_karyoplotter.py
# Author:     Marco A. Coelho @ Heitman lab, Duke University
# Date:       2025-04-16
# Version:    1.1
# Description:
#     Perform nucleotide-based synteny analysis between genome assemblies.
#     Runs pairwise comparisons using minimap2 (asm20 preset), parses synteny blocks,
#     and generates high-quality SVG plots. Optionally accepts a fixed reference genome
#     to generate a single composite figure with all query genomes aligned to that reference.
#
#     Output visualizations display chromosome-scale alignments, color-coded synteny segments,
#     and optionally annotated centromere positions. When available, centromeres are highlighted
#     as black bars with overlaid white ellipses to facilitate the detection of inter-centromeric
#     rearrangements and structural variation.
#
#     This version uses the minimap2 'asm20' preset, which is designed for genome pairs
#     with up to ~10% sequence divergence. In practice, we found it performs robustly
#     even for comparisons with divergence exceeding this threshold.
#
# Requirements:
#     - Python 3
#     - minimap2 (in $PATH, tested with version 2.17-r941 )
#     - Biopython
#     - pandas
#     - matplotlib
#     - numpy
#
# Input formats:
#     Genome assemblies must be in FASTA format (.fa or .fasta).
#     Each contig or chromosome must have a unique identifier in the FASTA header:
#         >chr_1
#         ATCGATCG...
#         >chr_2
#         GCTAGCTA...
#
#     Optional centromere coordinates should be provided as BED files with the same
#     prefix as the genome file:
#           Cryptococcus_neoformans_125.91.scaffolds.fa  →  Cryptococcus_neoformans_125.91.cen.bed
#
#     BED format:
#           chr_1   960921  1030510 CEN1
#           chr_2   904508  935491  CEN2
#           ...
#
# Usage:
#     python3 0_synteny_karyoplotter.py --input_dir /path/to/genomes
#     Optional: --ref genome_prefix (to compare all other genomes against this one)
# ==============================================================================



import os
import re
import argparse
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
from collections import defaultdict
from matplotlib.patches import Ellipse
from Bio import SeqIO
import numpy as np
from collections import OrderedDict


plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'


def parse_args():
    parser = argparse.ArgumentParser(description="Synteny analysis: run nucleotide comparisons and plot to-scale alignments.")
    parser.add_argument('--input_dir', required=True, help="Folder with genome files (.fa and optional .bed for centromeres)")
    parser.add_argument('--ref', help="Optional reference genome name (to compare against all others)")
    return parser.parse_args()

def get_max_contig_length_across_all_fastas(input_dir):
    """Scan all FASTA files and return the maximum contig length found."""
    max_length = 0
    for file in os.listdir(input_dir):
        if file.endswith((".fa", ".fasta")):
            for record in SeqIO.parse(os.path.join(input_dir, file), "fasta"):
                max_length = max(max_length, len(record.seq))
    return max_length


def get_contig_lengths(fasta_file):
    """Return a dict {contig: length} preserving FASTA order."""
    contig_lens = OrderedDict()
    for record in SeqIO.parse(fasta_file, "fasta"):
        contig_lens[record.id] = len(record.seq)
    return contig_lens


def normalize_name(filename):
    """Strip common suffixes and extensions to unify genome names."""
    name = filename
    name = re.sub(r'\.(fa|fasta|gff3|gbk|gbff)$', '', name)
    name = re.sub(r'\.(scaffolds|genome|assembly)$', '', name)
    return name

def scan_genomes(input_dir):
    genome_dict = defaultdict(dict)
    for file in os.listdir(input_dir):
        full_path = os.path.join(input_dir, file)
        if file.endswith((".fa", ".fasta", ".gff3", ".gbk", ".gbff")):
            name = normalize_name(file)
            if file.endswith((".fa", ".fasta")):
                genome_dict[name]['fasta'] = full_path
            elif file.endswith(".gff3"):
                genome_dict[name]['gff3'] = full_path
            elif file.endswith((".gbk", ".gbff")):
                genome_dict[name]['gbk'] = full_path
    return genome_dict

def generate_pairs(genome_dict, ref=None):
    names = list(genome_dict.keys())
    if ref and ref in names:
        return [(ref, other) for other in names if other != ref]
    else:
        return [(a, b) for a in names for b in names if a != b]

def run_minimap2(ref_fasta, qry_fasta, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    ref_base = normalize_name(os.path.basename(ref_fasta).rsplit('.', 1)[0])
    qry_base = normalize_name(os.path.basename(qry_fasta).rsplit('.', 1)[0])
    out_file = os.path.join(output_dir, f"{qry_base}_vs_{ref_base}.paf")
    
    cmd = ["minimap2", "-x", "asm20", ref_fasta, qry_fasta]
    with open(out_file, "w") as out:
        print(f"🔄 Running minimap2: {qry_base} → {ref_base}")
        subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE)
    
    print(f"✅ Output saved to: {out_file}")
    return out_file

def parse_paf(paf_file, output_dir="parsed_blocks"):
    os.makedirs(output_dir, exist_ok=True)
    columns = [
        "query", "query_len", "query_start", "query_end",
        "strand", "target", "target_len", "target_start", "target_end",
        "match_len", "block_len", "mapq"
    ]
    df = pd.read_csv(paf_file, sep="\t", header=None, usecols=range(12), names=columns)
    df["block_id"] = [f"block_{i+1}" for i in range(len(df))]

    out_path = os.path.join(output_dir, os.path.basename(paf_file).replace(".paf", "_synteny.tsv"))
    df.to_csv(out_path, sep="\t", index=False)
    print(f"🧬 Parsed synteny blocks saved to: {out_path}")
    return out_path

def read_centromeres(bed_file):
    centromeres = {}
    if not os.path.exists(bed_file):
        return centromeres  # Return empty if file doesn't exist

    with open(bed_file) as f:
        for line in f:
            if line.strip() == "" or line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                chrom, start, end = parts[:3]
                try:
                    centromeres[chrom] = (int(start), int(end))
                except ValueError:
                    continue  # Skip malformed lines
    return centromeres


def plot_synteny_blocks_final(tsv_file, ref_fasta, qry_fasta, ref_cen_bed=None, qry_cen_bed=None, output_dir="plots", global_max=None):
    os.makedirs(output_dir, exist_ok=True)

    df = pd.read_csv(tsv_file, sep="\t")
    ref_lengths = get_contig_lengths(ref_fasta)
    qry_lengths = get_contig_lengths(qry_fasta)

    ref_cens = {}
    qry_cens = {}

    if ref_cen_bed and os.path.exists(ref_cen_bed):
        ref_cens = read_centromeres(ref_cen_bed)
        print(f"📍 Loaded {len(ref_cens)} centromeres from {ref_cen_bed}")
    else:
        print(f"⚠️  No centromeres found for REF genome: {os.path.basename(ref_fasta)}")

    if qry_cen_bed and os.path.exists(qry_cen_bed):
        qry_cens = read_centromeres(qry_cen_bed)
        print(f"📍 Loaded {len(qry_cens)} centromeres from {qry_cen_bed}")
    else:
        print(f"⚠️  No centromeres found for QRY genome: {os.path.basename(qry_fasta)}")

    # Sort chromosomes
    ref_contigs = list(ref_lengths.keys())
    qry_contigs = list(qry_lengths.keys())

    max_length = global_max if global_max else max(max(ref_lengths.values()), max(qry_lengths.values()))
    buffer = int(max_length * 0.05)
    max_plot = max_length + buffer

    unit = 1_000_000 if max_length > 1_000_000 else 1_000
    ylabel = "Length (Mb)" if unit == 1_000_000 else "Length (kb)"

    # Use tab20 for color clarity
    cmap = plt.colormaps["tab20"]
    colors = [cmap(i % cmap.N) for i in range(len(ref_contigs))]
    color_map = {contig: colors[i] for i, contig in enumerate(ref_contigs)}


    # Custom fixed palette for first 14 chromosomes (you can modify these hex codes!)
    fixed_colors = [
        "#c4a5ce", "#f9df00", "#ec9bc2", "#f48130", "#dfb08c",
        "#e11f26", "#6b3e97", "#1a75bb", "#b31f5d", "#92d7f0",
        "#a7d069", "#31b7af", "#b05a28", "#34a148"
    ]

    # Fallback colormap for chromosomes beyond 14
    fallback_cmap = plt.get_cmap("tab20")

    colors = []
    for i in range(len(ref_contigs)):
        if i < len(fixed_colors):
            colors.append(fixed_colors[i])
        else:
            colors.append(fallback_cmap((i - len(fixed_colors)) % fallback_cmap.N))

    color_map = {contig: colors[i] for i, contig in enumerate(ref_contigs)}

    fig, ax = plt.subplots(figsize=(10, 8))
    bar_width = 0.4
    spacing_within = 0.7     # tighter spacing *within* each genome
    spacing_between = 2      # space *between* reference and query sides

    fixed_ellipse_height = 80000  # Set centromere ellipse symbol size here, used in both ref and query loops

    # Plot reference chromosomes (left)
    for i, contig in enumerate(ref_contigs):
        x = i * spacing_within
        height = ref_lengths[contig]
        ax.add_patch(patches.Rectangle((x - bar_width/2, 0), bar_width, height,
                                       facecolor=color_map[contig], edgecolor='black'))
        
        ax.text(x, -buffer * 0.6, str(i+1), ha='center', va='top', fontsize=14)

        if contig in ref_cens:
            start, end = ref_cens[contig]
            # Black centromere bar
            ax.add_patch(patches.Rectangle((x - bar_width/2, start), bar_width, end - start,
                                        facecolor='black', alpha=0.85, edgecolor='none', zorder=10))
            # Add white ellipse as centromere symbol (fixed size)
            cent_mid = start + (end - start) / 2
            ax.add_patch(Ellipse((x, cent_mid),
                                width=bar_width * 0.6,
                                height=fixed_ellipse_height,
                                facecolor='white', edgecolor='black',
                                lw=1.5, zorder=11))
            

    # Plot query chromosomes (right)
    x_offset = spacing_within * (len(ref_contigs) - 1) + spacing_between
    for i, contig in enumerate(qry_contigs):
        x = x_offset + i * spacing_within
        height = qry_lengths[contig]

        # Draw chromosome label
        ax.text(x, -buffer * 0.6, str(i+1), ha='center', va='top', fontsize=14)

        if contig in qry_cens:
            start, end = qry_cens[contig]
            # Black centromere bar (underneath the ellipse)
            ax.add_patch(patches.Rectangle((x - bar_width/2, start), bar_width, end - start,
                                        facecolor='black', alpha=0.85, edgecolor='none', zorder=10))

            # White ellipse symbol (centromere)
            cent_mid = start + (end - start) / 2
            ax.add_patch(Ellipse((x, cent_mid),
                                width=bar_width * 0.6,
                                height=fixed_ellipse_height,
                                facecolor='white', edgecolor='black',
                                lw=1.5, zorder=11))

        # Draw query chromosome border on top (no fill)
        ax.add_patch(patches.Rectangle((x - bar_width/2, 0), bar_width, height,
                                    facecolor='none', edgecolor='black', lw=1.2, zorder=20))


    # Paint query synteny segments based on reference color
    for _, row in df.iterrows():
        q_contig = row["query"]
        t_contig = row["target"]
        q_start = min(row["query_start"], row["query_end"])
        q_end = max(row["query_start"], row["query_end"])
        color = color_map.get(t_contig, "gray")

        if q_contig in qry_contigs:
            i = qry_contigs.index(q_contig)
            x = x_offset + i * spacing_within
            ax.add_patch(patches.Rectangle((x - bar_width/2, q_start), bar_width, q_end - q_start,
                                           facecolor=color, edgecolor=None, alpha=0.95))

    # Styling
    ax.invert_yaxis()
    ax.set_xlim(-1, x_offset + len(qry_contigs) * spacing_within + 1)
    ax.set_ylim(max_plot, -buffer)
    ax.set_ylabel(ylabel, fontsize=16)  # You can increase 14 to 16 or more


    ref_base = os.path.basename(ref_fasta).replace(".scaffolds.fa", "").replace(".fasta", "")
    qry_base = os.path.basename(qry_fasta).replace(".scaffolds.fa", "").replace(".fasta", "")

    # Label each species above its chromosome block
    mid_ref = (len(ref_contigs) - 1) * spacing_within / 2
    mid_qry = x_offset + (len(qry_contigs) - 1) * spacing_within / 2

    ax.text(mid_ref, -buffer * 1.8, ref_base, ha='center', va='top', fontsize=14, fontweight='regular')
    ax.text(mid_qry, -buffer * 1.8, qry_base, ha='center', va='top', fontsize=14, fontweight='regular')

    ax.set_xticks([])
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{x/unit:.1f}"))
    ax.tick_params(axis='y', labelsize=14)  # Adjust number size (default ~10)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.yaxis.grid(True, linestyle='--', alpha=0.3)

    # Save as SVG
    svg_file = os.path.join(output_dir, f"{ref_base}_vs_{qry_base}_final.svg")
    plt.tight_layout()
    plt.savefig(svg_file)
    print(f"🖼️ Final synteny figure saved to: {svg_file}")


def plot_synteny_multiple_queries(ref_fasta, qry_data, ref_cen_bed=None, global_max=None, output_dir="plots"):
    def get_contig_lengths(fasta_file):
        return {rec.id: len(rec.seq) for rec in SeqIO.parse(fasta_file, "fasta")}

    def read_centromeres(bed_file):
        cens = {}
        if not os.path.exists(bed_file):
            return cens
        with open(bed_file) as f:
            for line in f:
                if line.strip() and not line.startswith("#"):
                    parts = line.strip().split("\t")
                    if len(parts) >= 3:
                        try:
                            cens[parts[0]] = (int(parts[1]), int(parts[2]))
                        except ValueError:
                            continue
        return cens

    os.makedirs(output_dir, exist_ok=True)

    ref_lengths = get_contig_lengths(ref_fasta)
    ref_contigs = list(ref_lengths.keys())
    ref_cens = read_centromeres(ref_cen_bed) if ref_cen_bed else {}

    max_length = global_max if global_max else max(ref_lengths.values())
    buffer = int(max_length * 0.05)
    max_plot = max_length + buffer

    unit = 1_000_000 if max_length > 1_000_000 else 1_000
    ylabel = "Length (Mb)" if unit == 1_000_000 else "Length (kb)"


    # Custom fixed palette for first 14 chromosomes (you can modify these hex codes!)
    fixed_colors = [
        "#c4a5ce", "#f9df00", "#ec9bc2", "#f48130", "#dfb08c",
        "#e11f26", "#6b3e97", "#1a75bb", "#b31f5d", "#92d7f0",
        "#a7d069", "#31b7af", "#b05a28", "#34a148"
    ]

    # Fallback colormap for chromosomes beyond 14
    fallback_cmap = plt.get_cmap("tab20")

    colors = []
    for i in range(len(ref_contigs)):
        if i < len(fixed_colors):
            colors.append(fixed_colors[i])
        else:
            colors.append(fallback_cmap((i - len(fixed_colors)) % fallback_cmap.N))

    color_map = {contig: colors[i] for i, contig in enumerate(ref_contigs)}

    bar_width = 0.4
    spacing_within = 0.7    # tighter spacing *within* each genome
    spacing_between = 2     # space *between* reference and query sides
    
    fixed_ellipse_height = 80000  # Set centromere ellipse symbol size here here, used in both ref and query loops

    # Count total chromosomes to estimate width
    total_chrs = len(ref_contigs) + sum(len(pd.read_csv(entry["tsv_file"], sep="\t")['query'].unique()) for entry in qry_data)
    width = 2 + total_chrs * spacing_within * 0.7  # adjust multiplier for padding

    # Create dynamic figure
    fig, ax = plt.subplots(figsize=(width, 8))
    x_ref = 0

    # Plot reference chromosomes (left side)
    for i, contig in enumerate(ref_contigs):
        x = x_ref + i * spacing_within
        height = ref_lengths[contig]
        ax.add_patch(patches.Rectangle((x - bar_width/2, 0), bar_width, height,
                                       facecolor=color_map[contig], edgecolor='black'))
        
        # Draw chromosome label
        ax.text(x, -buffer * 0.6, str(i+1), ha='center', va='top', fontsize=14)

        if contig in ref_cens:
            start, end = ref_cens[contig]
            # Black centromere bar
            ax.add_patch(patches.Rectangle((x - bar_width/2, start), bar_width, end - start,
                                           facecolor='black', alpha=0.85, edgecolor='none', zorder=10))
            # Add white ellipse as centromere symbol (fixed size)
            cent_mid = start + (end - start) / 2
            ax.add_patch(Ellipse((x, cent_mid),
                                width=bar_width * 0.6,
                                height=fixed_ellipse_height,
                                facecolor='white', edgecolor='black',
                                lw=1.5, zorder=11))
        
    x_cursor = x_ref + (len(ref_contigs) - 1) * spacing_within + spacing_between

    # Plot all the query chromosomes (right side)
    for entry in qry_data:
        df = pd.read_csv(entry["tsv_file"], sep="\t")
        qry_lengths = get_contig_lengths(entry["qry_fasta"])
        qry_contigs = list(qry_lengths.keys())
        qry_cens = read_centromeres(entry["qry_cen_bed"]) if os.path.exists(entry["qry_cen_bed"]) else {}

        for i, contig in enumerate(qry_contigs):
            x = x_cursor + i * spacing_within
            height = qry_lengths[contig]
            ax.add_patch(patches.Rectangle((x - bar_width/2, 0), bar_width, height,
                                           facecolor='none', edgecolor='none', linewidth=1.2))
            
            # Draw chromosome label
            ax.text(x, -buffer * 0.6, str(i+1), ha='center', va='top', fontsize=14)

            if contig in qry_cens:
                start, end = qry_cens[contig]
                # Black centromere bar (underneath the ellipse)
                ax.add_patch(patches.Rectangle((x - bar_width/2, start), bar_width, end - start,
                                               facecolor='black', alpha=0.85, edgecolor='none', zorder=10))
                # White ellipse symbol (centromere)
                cent_mid = start + (end - start) / 2
                ax.add_patch(Ellipse((x, cent_mid),
                                width=bar_width * 0.6,
                                height=fixed_ellipse_height,
                                facecolor='white', edgecolor='black',
                                lw=1.5, zorder=11))    
            
            # Draw query chromosome border on top (no fill)
            ax.add_patch(patches.Rectangle((x - bar_width/2, 0), bar_width, height,
                                    facecolor='none', edgecolor='black', lw=1.2, zorder=20))

        for _, row in df.iterrows():
            q_contig = row["query"]
            t_contig = row["target"]
            q_start = min(row["query_start"], row["query_end"])
            q_end = max(row["query_start"], row["query_end"])
            color = color_map.get(t_contig, "gray")
            if q_contig in qry_contigs:
                i = qry_contigs.index(q_contig)
                x = x_cursor + i * spacing_within
                ax.add_patch(patches.Rectangle((x - bar_width/2, q_start), bar_width, q_end - q_start,
                                               facecolor=color, edgecolor=None, alpha=0.95))

        # Label species name
        qry_label = os.path.basename(entry["qry_fasta"]).replace(".scaffolds.fa", "").replace(".fasta", "")
        mid_x = x_cursor + (len(qry_contigs) - 1) * spacing_within / 2
        ax.text(mid_x, -buffer * 1.8, qry_label, ha='center', va='top', fontsize=14, fontweight='regular')

        x_cursor += len(qry_contigs) * spacing_within + spacing_between

    # Final layout
    ref_base = os.path.basename(ref_fasta).replace(".scaffolds.fa", "").replace(".fasta", "")
    mid_ref = (len(ref_contigs) - 1) * spacing_within / 2
    ax.text(mid_ref, -buffer * 1.8, ref_base, ha='center', va='top', fontsize=14, fontweight='regular')

    ax.invert_yaxis()
    ax.set_xlim(-1, x_cursor)
    ax.set_ylim(max_plot, -buffer)
    ax.set_ylabel(ylabel, fontsize=16)      # Y-axis font label
    ax.set_xticks([])
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{x/unit:.1f}"))
    ax.tick_params(axis='y', labelsize=14)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.yaxis.grid(True, linestyle='--', alpha=0.3)

    plt.tight_layout()
    out_path = os.path.join(output_dir, f"{ref_base}_multi_query_synteny.svg")
    try:
        plt.tight_layout()
        plt.savefig(out_path)
        print(f"🖼️ Multi-query synteny figure saved to: {out_path}")
    except Exception as e:
        print(f"❌ Failed to save multi-query plot: {e}")


def main():
    args = parse_args()
    genomes = scan_genomes(args.input_dir)
    global_max_length = get_max_contig_length_across_all_fastas(args.input_dir)
    print(f"📏 Global max contig length across all genomes: {global_max_length}")

    print(f"🔍 Found {len(genomes)} genome(s) in '{args.input_dir}'")
    for name, files in genomes.items():
        print(f"  📦 {name}: {list(files.keys())}")

    if args.ref:
        ref = args.ref
        if ref not in genomes:
            print(f"❌ Reference genome '{ref}' not found in input folder.")
            return

        print(f"📌 Using fixed reference genome: {ref}")
        ref_info = genomes[ref]
        ref_base = normalize_name(os.path.basename(ref_info['fasta']).replace(".scaffolds", ""))
        ref_cen_file = os.path.join(args.input_dir, f"{ref_base}.cen.bed")

        queries = [g for g in genomes if g != ref]
        all_data = []

        print(f"\n📊 Preparing {len(queries)} comparisons vs reference: {ref}")

        for qry in queries:
            qry_info = genomes[qry]
            if 'fasta' not in qry_info:
                print(f"⚠️ Missing FASTA for query genome: {qry}")
                continue

            qry_base = normalize_name(os.path.basename(qry_info['fasta']).replace(".scaffolds", ""))
            qry_cen_file = os.path.join(args.input_dir, f"{qry_base}.cen.bed")

            print(f"✅ Nucleotide: {qry} (query) vs {ref} (reference)")
            paf_file = run_minimap2(ref_info['fasta'], qry_info['fasta'], output_dir="results")
            tsv_file = parse_paf(paf_file)

            all_data.append({
                "tsv_file": tsv_file,
                "qry_fasta": qry_info['fasta'],
                "qry_cen_bed": qry_cen_file,
                "qry_name": qry
            })

        print(f"🧪 Plotting {len(all_data)} genome(s) against reference '{ref}'")
        plot_synteny_multiple_queries(
            ref_fasta=ref_info['fasta'],
            qry_data=all_data,
            ref_cen_bed=ref_cen_file,
            global_max=global_max_length,
            output_dir="plots"
        )

    else:
        pairs = generate_pairs(genomes)
        print(f"\n📊 Preparing {len(pairs)} pairwise comparisons:")
        if not pairs:
            print("⚠️  No valid genome pairs found. Check FASTA file naming or normalization.")
            return

        for ref, qry in pairs:
            ref_info = genomes[ref]
            qry_info = genomes[qry]

            if 'fasta' in ref_info and 'fasta' in qry_info:
                print(f"✅ Nucleotide: {qry} (query) vs {ref} (reference)")
                paf_file = run_minimap2(ref_info['fasta'], qry_info['fasta'], output_dir="results")
                tsv_file = parse_paf(paf_file)

                ref_base = normalize_name(os.path.basename(ref_info['fasta']).replace(".scaffolds", ""))
                qry_base = normalize_name(os.path.basename(qry_info['fasta']).replace(".scaffolds", ""))
                ref_cen_file = os.path.join(args.input_dir, f"{ref_base}.cen.bed")
                qry_cen_file = os.path.join(args.input_dir, f"{qry_base}.cen.bed")

                print(f"🔍 Looking for centromere files:\n   REF: {ref_cen_file}\n   QRY: {qry_cen_file}")

                plot_synteny_blocks_final(
                    tsv_file,
                    ref_info['fasta'],
                    qry_info['fasta'],
                    ref_cen_bed=ref_cen_file,
                    qry_cen_bed=qry_cen_file,
                    global_max=global_max_length
                )
            else:
                print(f"❌ Missing FASTA for pair: {qry} vs {ref}")

    print("🧩 All comparisons completed.")


if __name__ == '__main__':
    main()