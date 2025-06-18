# 0_plot_dDDH_outputs_from_GGDC.py

# ==============================================================================
# Author: Marco A. Coelho @ Heitman Lab, Duke University
# Date: 2025-04-14
# Version: 1.9
# Description: 
#     Parses GGDC results (CSV with headers on second row), extracts Formula 2 DDH values,
#     classifies intra/inter-species comparisons, generates clean strain labels,
#     and produces selected plots and summary outputs based on user-defined subset.
# Requirements:
#     - pandas
#     - matplotlib
#     - seaborn
# Usage:
#     python3 0_plot_dDDH_outputs_from_GGDC.py <ggdc_csv_input>
# ==============================================================================

import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch

# --- Global settings ---
plt.rcParams['pdf.fonttype'] = 42
sns.set_theme(style="whitegrid")

# --- Font size settings ---
title_fontsize = 14
label_fontsize = 12
tick_fontsize = 10
legend_fontsize = 10

# --- Classification function ---
def classify_species_relationship(ddh, prob):
    if ddh >= 70:
        return "Same species (by DDH)"
    elif ddh >= 67:
        if prob >= 90:
            return "Likely same species (high confidence)"
        elif prob >= 70:
            return "Likely same species (moderate)"
        elif prob >= 50:
            return "Borderline / uncertain"
    return "Different species"

# --- Argument parsing ---
if len(sys.argv) != 2:
    print("Usage: python3 0_plot_dDDH_outputs_from_GGDC.py <ggdc_csv_input>")
    sys.exit(1)

input_csv = sys.argv[1]
if not os.path.exists(input_csv):
    print(f"❌ Input file not found: {input_csv}")
    sys.exit(1)

# --- Derive file prefix ---
prefix = os.path.splitext(os.path.basename(input_csv))[0]

# --- Load data ---
df = pd.read_csv(input_csv, header=1)
df["DDH_F2"] = pd.to_numeric(df["DDH.1"], errors="coerce")
df["Prob_ge_70"] = pd.to_numeric(df["Prob. DDH >= 70%.1"], errors="coerce")

df["Species_1"] = df["Query genome"].apply(lambda x: "_".join(os.path.basename(x).split("_")[0:2]))
df["Species_2"] = df["Reference genome"].apply(lambda x: "_".join(os.path.basename(x).split("_")[0:2]))
df["Type"] = df.apply(lambda row: "Intra" if row["Species_1"] == row["Species_2"] else "Inter", axis=1)

def clean_strain_label(row):
    def clean(name):
        base = os.path.basename(name).replace(".scaffolds.fa", "").replace(".ENR", "")
        parts = base.split("_")
        return parts[-1] if len(parts) > 1 else base
    return f"{clean(row['Query genome'])} vs {clean(row['Reference genome'])}"

df["CleanLabel"] = df.apply(clean_strain_label, axis=1)

summary = df[["CleanLabel", "Species_1", "Species_2", "Type", "DDH_F2", "Prob_ge_70"]].dropna()
summary = summary.sort_values(by="DDH_F2", ascending=False)
summary["Classification"] = summary.apply(lambda row: classify_species_relationship(row["DDH_F2"], row["Prob_ge_70"]), axis=1)

summary_file = f"{prefix}_formula2_summary_strain_labels.csv"
summary.to_csv(summary_file, index=False)
print(f"✔ Saved summary table: {summary_file}")

color_map = {"Intra": "#ef771d", "Inter": "#1f77b4"}

# === BAR PLOT ===
plt.figure(figsize=(10, 6))
ax = sns.barplot(data=summary, x="CleanLabel", y="DDH_F2", hue="Type", palette=color_map, dodge=False, edgecolor="none")
plt.axhline(70, color='grey', linestyle='--')
ax.set_title("Genome-to-genome distance (formula 2) by strain comparison", fontsize=title_fontsize)
ax.set_ylabel("Formula 2 DDH (%)", fontsize=label_fontsize)
ax.set_xlabel("")
ax.tick_params(axis='x', labelsize=tick_fontsize, rotation=90)
ax.legend(handles=[
    Patch(facecolor="#ef771d", label="Intra-species"),
    Patch(facecolor="#1f77b4", label="Inter-species"),
    plt.Line2D([0], [0], color='gray', linestyle='--', label="70% threshold")
], fontsize=legend_fontsize)
plt.tight_layout()
plt.savefig(f"{prefix}_ddh_formula2_barplot_strain_labels.pdf", dpi=300)
plt.close()

# === SCATTER PLOT: DDH vs probability ===
plt.figure(figsize=(5, 5))
ax = sns.scatterplot(data=summary, x="DDH_F2", y="Prob_ge_70", hue="Type", palette=color_map, s=50, edgecolor="none", alpha=0.7)
plt.axvline(70, ls='--', color='gray')
plt.axhline(70, ls='--', color='gray')
ax.set_xlabel("Formula 2 DDH (%)", fontsize=label_fontsize)
ax.set_ylabel("Probability DDH ≥ 70%", fontsize=label_fontsize)
ax.set_title("DDH vs. probability of species identity", fontsize=title_fontsize)
plt.tight_layout()
plt.savefig(f"{prefix}_ddh_vs_probability_scatter_strain_labels.pdf", dpi=300)
plt.close()

# === SCATTER PLOT: intra-species, colored by species ===
intra = summary[summary["Type"] == "Intra"].copy()
plt.figure(figsize=(6, 6))
ax = sns.scatterplot(data=intra, x="DDH_F2", y="Prob_ge_70", hue="Species_1", s=50, edgecolor="none", alpha=0.7)
plt.axvline(70, ls='--', color='gray')
plt.axhline(70, ls='--', color='gray')
ax.set_title("DDH vs. probability (intra-species, colored by species)", fontsize=title_fontsize)
ax.set_xlabel("Formula 2 DDH (%)", fontsize=label_fontsize)
ax.set_ylabel("Probability DDH ≥ 70%", fontsize=label_fontsize)
plt.tight_layout()
plt.savefig(f"{prefix}_scatter_intra_by_species.pdf", dpi=300)
plt.close()

# === SCATTER PLOT: inter-species comparisons with labels ===
inter = summary[summary["Type"] == "Inter"].copy().reset_index(drop=True)
inter["Pair_ID"] = inter.index + 1

plt.figure(figsize=(6, 6))
ax = sns.scatterplot(data=inter, x="DDH_F2", y="Prob_ge_70", color="#1f77b4", s=60, edgecolor="none", alpha=0.7)
for i, row in inter.iterrows():
    ax.text(row["DDH_F2"] + 0.5, row["Prob_ge_70"], str(row["Pair_ID"]), fontsize=8, color='black')
plt.axvline(70, ls='--', color='gray')
plt.axhline(70, ls='--', color='gray')
ax.set_title("DDH vs. probability (inter-species comparisons)", fontsize=title_fontsize)
ax.set_xlabel("Formula 2 DDH (%)", fontsize=label_fontsize)
ax.set_ylabel("Probability DDH ≥ 70%", fontsize=label_fontsize)
plt.tight_layout()
plt.savefig(f"{prefix}_scatter_inter_numbered.pdf", dpi=300)
plt.close()

# === Export inter-species pair legend ===
inter_legend = inter[["Pair_ID", "Species_1", "Species_2", "CleanLabel", "DDH_F2", "Prob_ge_70"]]
inter_legend_file = f"{prefix}_inter_species_pair_legend.csv"
inter_legend.to_csv(inter_legend_file, index=False)
print(f"✔ Saved inter-species pair legend table: {inter_legend_file}")
