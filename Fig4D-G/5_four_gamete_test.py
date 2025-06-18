#!/usr/bin/env python3

#===============================================================================
# Author: Marco A. Coelho @ Heitman Lab, Duke University
# Date: 2025-05-29
# Version: 1.0
#
# Description:
#   Performs a four-gamete test on biallelic SNPs from a multi-sample VCF file
#   from Snippy. Measures the proportion of SNP pairs that violate
#   the four-gamete condition across distance bins — a signal of historical
#   recombination.
#
# Requirements:
#   - Conda environment: popgen_env
#   - Python libraries: pyvcf, pandas, matplotlib, tqdm
#
# Usage:
#   conda activate popgen_env
#   python 5_four_gamete_test.py
#
# Input:
#   core.vcf                        # Multi-sample, biallelic SNP VCF
#
# Output:
#   four_gamete_test_plot.pdf      # Barplot summarizing violations by distance bin
#   four_gamete_test_summary.tsv   # TSV table with number and percent of violating SNP pairs per bin
#===============================================================================

import vcf
from collections import defaultdict
import pandas as pd
import matplotlib.pyplot as plt
from itertools import combinations
from tqdm import tqdm

# --- PLOTTING SETTINGS ---
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'

# --- PARAMETERS ---
vcf_path = "core.vcf"
max_distance = 10000
distance_bins = [(0, 1000), (1000, 5000), (5000, 10000)]

# --- PARSE VCF AND EXTRACT BIALLELIC SNPs ---
vcf_reader = vcf.Reader(filename=vcf_path)
snps = []

for record in vcf_reader:
    if record.is_snp and len(record.alleles) == 2:
        haplotypes = [call.gt_bases for call in record.samples]
        if None in haplotypes or any(gt is None for gt in haplotypes):
            continue

        alleles = []
        for call in record.samples:
            gt = call.gt_bases
            if gt is None or call.is_het:
                break
            allele = gt.replace(record.REF, "0").replace(str(record.ALT[0]), "1")
            alleles.append(allele)
        if len(alleles) != len(record.samples):
            continue

        snps.append((record.CHROM, record.POS, alleles))

# --- FOUR-GAMETE TEST ---
violations_by_bin = defaultdict(int)
total_by_bin = defaultdict(int)

for (i, (chrom1, pos1, hap1)), (j, (chrom2, pos2, hap2)) in tqdm(
    combinations(enumerate(snps), 2),
    total=len(snps) * (len(snps) - 1) // 2,
    desc="Running four-gamete test"
):
    if chrom1 != chrom2:
        continue
    distance = abs(pos1 - pos2)
    if distance > max_distance:
        continue

    hap_pairs = set(zip(hap1, hap2))
    for bin_start, bin_end in distance_bins:
        if bin_start <= distance < bin_end:
            total_by_bin[(bin_start, bin_end)] += 1
            if len(hap_pairs) == 4:
                violations_by_bin[(bin_start, bin_end)] += 1

# --- SUMMARIZE RESULTS ---
results = []
for bin_range in distance_bins:
    total = total_by_bin[bin_range]
    violations = violations_by_bin[bin_range]
    percent = 100 * violations / total if total > 0 else 0
    results.append({
        "Distance bin": f"{bin_range[0]}–{bin_range[1]} bp",
        "SNP pairs": total,
        "Violations": violations,
        "% Violating": round(percent, 2)
    })

df_results = pd.DataFrame(results)
print("\n Four-gamete test summary:\n")
print(df_results.to_string(index=False))

# Save table to TSV
df_results.to_csv("four_gamete_test_summary.tsv", sep="\t", index=False)

# --- PLOT RESULTS ---
plt.figure(figsize=(6, 4))
plt.bar(
    df_results["Distance bin"],
    df_results["% Violating"],
    color="skyblue",
    edgecolor="black"
)
plt.ylabel("SNP pairs violating 4-gamete test (%)")
plt.xlabel("Distance between SNP pairs")
plt.title("Evidence of recombination (four-gamete test)")
plt.tight_layout()
plt.savefig("four_gamete_test_plot.pdf")
plt.show()
