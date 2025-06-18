#!/usr/bin/env python3

#===============================================================================
# Author: Marco A. Coelho @ Heitman Lab, Duke University
# Date: 2025-05-29
# Version: 1.0
#
# Description:
#   Visualizes LD decay from PopLDdecay output files with and without MAF filtering.
#   Generates smoothed r² plots (rolling average + LOWESS), calculates LD 50% decay 
#   distances, and exports key data tables and figure files.
#
# Requirements:
#   - Conda environment: popgen_env
#   - Python libraries: pandas, matplotlib, seaborn, statsmodels
#
# Usage:
#   conda activate popgen_env
#   python 4_plot_ld_decay.py
#
# Input:
#   core_lddecay_maf0.167_dist50000.stat.gz
#   core_lddecay_nomaf_dist50000.stat.gz
#
# Output:
#   ld_decay_maf.pdf / .png
#   ld_decay_nomaf.pdf / .png
#   ld_decay_combined.pdf / .png
#   ld_decay_lowess_maf.tsv
#   ld_decay_lowess_nomaf.tsv
#   ld_decay_summary.txt
#===============================================================================

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import gzip
import sys
from statsmodels.nonparametric.smoothers_lowess import lowess

# --- PLOTTING SETTINGS ---
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'
sns.set(style="whitegrid")

# --- LOGGING TO FILE AND STDOUT ---
class Tee:
    def __init__(self, *files):
        self.files = files
    def write(self, obj):
        for f in self.files:
            f.write(obj)
    def flush(self):
        for f in self.files:
            f.flush()

logfile = open("ld_decay_summary.txt", "w")
sys.stdout = Tee(sys.stdout, logfile)

# --- LOAD LD DECAY DATA ---
def load_lddecay(file_path, label):
    with gzip.open(file_path, 'rt') as f:
        df = pd.read_csv(f, sep='\t', comment='#')
        f.seek(0)
        header = f.readline().strip().lstrip('#').split('\t')
        df.columns = header
    df = df.rename(columns={'#Dist': 'Distance', 'Dist': 'Distance', 'Mean_r^2': 'R2'})
    df['Distance'] = pd.to_numeric(df['Distance'], errors='coerce')
    df['R2'] = pd.to_numeric(df['R2'], errors='coerce')
    df['Label'] = label
    return df[df['Distance'] <= 250000].sort_values("Distance")

df_maf = load_lddecay("core_lddecay_maf0.167_dist50000.stat.gz", "MAF ≥ 0.167")
df_nomaf = load_lddecay("core_lddecay_nomaf_dist50000.stat.gz", "No MAF filter")

# --- ROLLING AVERAGE SMOOTHING ---
df_maf["R2_smooth"] = df_maf["R2"].rolling(window=5, center=True, min_periods=1).mean()
df_nomaf["R2_smooth"] = df_nomaf["R2"].rolling(window=5, center=True, min_periods=1).mean()

# --- LOWESS SMOOTHING ---
df_maf_sub = df_maf.iloc[::10]
df_nomaf_sub = df_nomaf.iloc[::10]
lowess_maf = lowess(df_maf_sub["R2_smooth"], df_maf_sub["Distance"], frac=0.1, return_sorted=True)
lowess_nomaf = lowess(df_nomaf_sub["R2_smooth"], df_nomaf_sub["Distance"], frac=0.1, return_sorted=True)

# --- INDIVIDUAL PLOTS ---
def plot_individual(df, label, color, filename, lowess_line):
    plt.figure(figsize=(8, 5))
    sns.lineplot(data=df, x="Distance", y="R2_smooth", color=color, linewidth=2, label=label)
    plt.plot(lowess_line[:, 0], lowess_line[:, 1], color="black", linestyle="--", label="LOWESS trend")
    plt.title(f"LD Decay (r² vs distance ≤ 250 kb)\n{label}")
    plt.xlabel("Distance between SNPs (bp)")
    plt.ylabel("Smoothed mean r²")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{filename}.pdf")
    plt.savefig(f"{filename}.png", dpi=300)
    plt.close()

plot_individual(df_maf, "MAF ≥ 0.167", "orange", "ld_decay_maf", lowess_maf)
plot_individual(df_nomaf, "No MAF filter", "steelblue", "ld_decay_nomaf", lowess_nomaf)

# --- COMBINED PLOT ---
plt.figure(figsize=(8, 5))
sns.lineplot(data=df_maf, x="Distance", y="R2_smooth", color="orange", label="MAF ≥ 0.167")
sns.lineplot(data=df_nomaf, x="Distance", y="R2_smooth", color="steelblue", label="No MAF filter")
plt.plot(lowess_maf[:, 0], lowess_maf[:, 1], color="darkorange", linestyle="--", label="MAF trend")
plt.plot(lowess_nomaf[:, 0], lowess_nomaf[:, 1], color="navy", linestyle="--", label="No MAF trend")
plt.title("LD Decay (r² vs distance ≤ 250 kb)")
plt.xlabel("Distance between SNPs (bp)")
plt.ylabel("Smoothed mean r²")
plt.legend()
plt.tight_layout()
plt.savefig("ld_decay_combined.pdf")
plt.savefig("ld_decay_combined.png", dpi=300)
plt.close()

# --- EXPORT LOWESS RESULTS ---
pd.DataFrame(lowess_maf, columns=["Distance", "R2"]).to_csv("ld_decay_lowess_maf.tsv", sep="\t", index=False)
pd.DataFrame(lowess_nomaf, columns=["Distance", "R2"]).to_csv("ld_decay_lowess_nomaf.tsv", sep="\t", index=False)

# --- LD 50% DECAY DISTANCE ---
print("\n LD 50% decay distances (based on smoothed r²):")
for df, label in [(df_maf, "MAF ≥ 0.167"), (df_nomaf, "No MAF filter")]:
    max_r2 = df["R2_smooth"].max()
    half_r2 = max_r2 * 0.5
    decay = df[df["R2_smooth"] <= half_r2]
    if not decay.empty:
        dist_50 = decay.iloc[0]["Distance"]
        print(f"  {label}: r² dropped to 50% at ~{int(dist_50):,} bp")
    else:
        print(f"  {label}: r² did not drop to 50% within 250 kb")

# --- R2 VALUES AT SELECTED DISTANCES ---
target_distances = [1000, 5000, 10000, 20000, 50000]
print("\n Smoothed r² values at target distances:")
for df, label in [(df_maf, "MAF ≥ 0.167"), (df_nomaf, "No MAF filter")]:
    print(f"\n  {label}:")
    for d in target_distances:
        nearest = df.iloc[(df["Distance"] - d).abs().argsort()[:1]]
        dist_val = int(nearest["Distance"].values[0])
        r2_val = round(nearest["R2_smooth"].values[0], 4)
        print(f"    ~{dist_val:,} bp: r² = {r2_val}")

# --- CLEANUP ---
sys.stdout = sys.__stdout__
logfile.close()
