#!/usr/bin/env python3

# ===============================================================================
# Author: Marco A. Coelho @ Heitman lab, Duke University
# Date: 2025-05-30
# Version: 2.0
# Filename: biolog_phenotype_analysis.py
#
# Description:
#   Full version of Biolog phenotype analysis across Cryptococcus clades.
#   Performs statistical tests, dimensionality reduction, clustering,
#   feature selection via Random Forest, and data visualization.
#   Includes heatmap, clustermap, and confusion matrix outputs.
#
# Requirements:
#   pandas, openpyxl, matplotlib, seaborn, scipy, scikit-learn, statsmodels
#
# Recommended usage:
#   Create a conda environment from the provided environment.yml file:
#     conda env create -f environment.yml
#     conda activate biolog_env
#
# Input format:
#   The input Excel file should be a binary phenotype matrix where:
#     - Rows represent phenotypic tests (e.g., substrates or conditions).
#     - Columns represent Cryptococcus strains (e.g., A_Cneo_H99).
#     - The first column must be named or treated as "Substrate" and contain test names.
#     - Remaining columns must contain phenotype calls as "positive", "negative", or blank.
#     - Cells will be automatically lowercased and parsed as binary (positive=1, negative=0).
#     - Strain names must include clade prefix (e.g., "A_", "B_", etc.) for downstream clade-based analysis.
#
# Usage:
#   python biolog_phenotype_analysis.py Cryptococcus_Biolog_data_parsed.xlsx
# ===============================================================================



# -------------------- Imports -------------------- #
# System and file handling
import sys
import os

# Data manipulation and I/O
import pandas as pd

# Plotting and visualization
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from matplotlib.patches import Ellipse
from matplotlib import colormaps

# Math and statistics
import numpy as np
from scipy.stats import fisher_exact, chi2_contingency, chi2
from scipy.cluster.hierarchy import linkage, dendrogram

# Machine learning
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, confusion_matrix, ConfusionMatrixDisplay

# Multiple testing correction
from statsmodels.stats.multitest import multipletests

# Text normalization
from unidecode import unidecode

# Matplotlib font compatibility for export
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'


# -------------------- Section 1: Load and clean data -------------------- #

def load_and_clean_data(filepath):
    """
    Load Biolog phenotype matrix from Excel and convert to binary format (1/0).
    - Assumes first column is 'Substrate' (test names)
    - Converts 'positive'/'negative' calls to 1/0
    - Returns a binary DataFrame: rows = phenotypes, columns = strains
    """
    # Read Excel file into DataFrame
    df_raw = pd.read_excel(filepath)

    # Rename the first column to 'Substrate' and set as index
    df_raw.rename(columns={df_raw.columns[0]: "Substrate"}, inplace=True)
    df_raw.set_index("Substrate", inplace=True)

    # Normalize entries: lowercase, remove leading/trailing whitespace
    df_clean = df_raw.astype(str).apply(lambda col: col.str.lower().str.strip())

    # Replace text with binary values: positive → 1, negative → 0
    df_bin = df_clean.replace({"positive": 1, "negative": 0})

    # Ensure integer type and fill missing values with 0
    df_bin = df_bin.infer_objects(copy=False).fillna(0).astype(int)

    return df_bin

# -------------------- Section 2 – Plate Annotation and Helpers -------------------- #

def greekify(name):
    """
    Replace common ASCII prefixes with Greek symbols for nicer labels.
    e.g., 'alpha-' becomes 'α-', 'Beta-' becomes 'Β-', etc.
    """
    return (name.replace("alpha-", "α-")
                .replace("beta-", "β-")
                .replace("gamma-", "γ-")
                .replace("Alpha-", "Α-")
                .replace("Beta-", "Β-")
                .replace("Gamma-", "Γ-"))

def normalize_test_name(name):
    """
    Helper function to standardize phenotype test names:
    - lowercase
    - stripped whitespace
    - ASCII-normalized (removes accents)
    """
    return unidecode(name.lower().strip())

def annotate_phenotype_with_plate_from_well_table(df, tsv_file):
    """
    Append Biolog plate origin to phenotype test names (e.g., '(GENIII)').

    Input:
    - df: phenotype DataFrame (rows = phenotypes)
    - tsv_file: tab-separated file with 'Test', 'YT_Well', 'FF_Well', 'GENIII_Well' columns

    Output:
    - Annotated DataFrame with updated row names (e.g., 'L-Proline (GENIII)')
    """

    df_annotated = df.copy()
    df_annotated.index.name = "Substrate"

    # Read plate mapping from TSV
    ref_df = pd.read_csv(tsv_file, sep="\t", encoding="ISO-8859-1")
    plate_columns = ["YT_Well", "FF_Well", "GENIII_Well"]
    plate_labels = {"YT_Well": "YT", "FF_Well": "FF", "GENIII_Well": "GENIII"}

    # Build a lookup: normalized test name → list of plate names
    plate_map = {}
    for _, row in ref_df.iterrows():
        test = normalize_test_name(row["Test"])
        plates = [plate_labels[col] for col in plate_columns if pd.notna(row[col])]
        if plates:
            plate_map.setdefault(test, set()).update(plates)

    # Update each phenotype label with its plate(s) if available
    new_index = []
    for test in df_annotated.index:
        norm = normalize_test_name(test)
        plates = sorted(plate_map.get(norm, []))
        new_index.append(f"{test} ({'/'.join(plates)})" if plates else test)

    df_annotated.index = new_index
    return df_annotated


# -------------------- Section 3 – Heatmap and Clustering Visualization -------------------- #

def plot_heatmap(df, output_file):
    """
    Plot a basic heatmap of phenotype presence/absence (binary values).
    Rows = phenotypes, columns = strains.
    """
    plt.figure(figsize=(14, 10))
    sns.heatmap(
        df,
        cmap="viridis",
        cbar_kws={'label': 'Presence (1) / Absence (0)'},
        linewidths=0.3
    )
    plt.title("Phenotypic Profile of All Strains")
    plt.xlabel("Strain")
    plt.ylabel("Phenotypic Test")
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

def plot_dendrogram(df, output_file):
    """
    Generate a hierarchical clustering dendrogram for strains
    using average linkage and Euclidean distance.
    """
    df_cluster = df.fillna(0).T  # Strains as rows
    linkage_matrix = linkage(df_cluster, method="average", metric="euclidean")

    plt.figure(figsize=(12, 6))
    dendrogram(
        linkage_matrix,
        labels=df_cluster.index,
        leaf_rotation=90
    )
    plt.title("Hierarchical Clustering of Strains")
    plt.xlabel("Strain")
    plt.ylabel("Distance")
    plt.tight_layout()
    plt.savefig(output_file)
    plt.savefig(output_file.replace(".pdf", ".png"), dpi=300)
    plt.close()

def plot_clustermap(df, output_file):
    """
    Generate a clustermap with:
    - hierarchical clustering
    - clade color annotations
    - binary colormap
    """
    df_cluster = df.fillna(0)
    strain_labels = df_cluster.columns
    clades = [s.split("_")[0] for s in strain_labels]

    # Create color map for clades
    clade_palette = sns.color_palette("Set2", len(set(clades)))
    clade_lut = dict(zip(sorted(set(clades)), clade_palette))
    col_colors = pd.Series(clades, index=strain_labels).map(clade_lut)

    # Binary color map for heatmap
    binary_cmap = mcolors.ListedColormap(["#eae7d9", "#726658"])

    # Create seaborn clustermap
    g = sns.clustermap(
        df_cluster,
        cmap=binary_cmap,
        col_colors=col_colors,
        figsize=(10, 30),
        cbar_kws={'label': 'Phenotype', 'ticks': [0, 1], 'shrink': 0.2},
        linewidths=0.2,
        xticklabels=True,
        yticklabels=True
    )

    # Minor axis cleanup
    g.ax_heatmap.tick_params(axis='y', left=False)
    g.fig.suptitle("Clustermap of Phenotypic Profiles with Clade Annotation", y=1.03)
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90)

    # Save outputs
    plt.tight_layout()
    g.savefig(output_file)
    g.savefig(output_file.replace(".pdf", ".png"), dpi=300)
    plt.close()

def plot_clustermap_versions(df, prefix="output"):
    """
    Generate two clustermaps:
    - Full version with all phenotypes
    - Filtered version with only variable phenotypes
    Also saves filtered matrix to CSV.
    """
    # Full clustermap
    plot_clustermap(df, output_file=f"{prefix}_clustermap_full.pdf")

    # Filter out non-informative phenotypes
    df_filtered = df.loc[df.nunique(axis=1) > 1]
    df_filtered.to_csv(f"{prefix}_filtered_binary_matrix.csv")

    # Only plot filtered clustermap if it has multiple phenotypes
    if df_filtered.shape[0] >= 2:
        plot_clustermap(df_filtered, output_file=f"{prefix}_clustermap_filtered.pdf")


# -------------------- Section 4 – PCA Analysis with Confidence Ellipses -------------------- #

def plot_pca(df, output_file):
    """
    Perform PCA on strain phenotypic profiles and visualize in 2D space.
    - Includes 95% confidence ellipses per clade.
    - Labels each point with strain name.
    """
    # Transpose: rows = strains, columns = phenotypes
    data_pca = df.fillna(0).T
    strain_labels = data_pca.index
    clades = [s.split("_")[0] for s in strain_labels]

    # Standardize features before PCA
    scaler = StandardScaler()
    data_scaled = scaler.fit_transform(data_pca)

    # Compute first two principal components
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(data_scaled)

    # Store results in DataFrame
    pca_df = pd.DataFrame(pca_result, columns=["PC1", "PC2"])
    pca_df["Clade"] = clades
    pca_df["Strain"] = strain_labels

    # Start plotting
    plt.figure(figsize=(10, 7))
    ax = sns.scatterplot(
        data=pca_df,
        x="PC1", y="PC2",
        hue="Clade",
        s=100, edgecolor="black"
    )

    # Annotate each point with strain name
    for _, row in pca_df.iterrows():
        ax.text(row["PC1"] + 0.2, row["PC2"], row["Strain"], fontsize=8)

    # Draw 95% confidence ellipses for each clade
    for clade, group in pca_df.groupby("Clade"):
        if len(group) < 3:
            continue  # Need ≥3 points for ellipse

        cov = np.cov(group[["PC1", "PC2"]].T)
        mean = group[["PC1", "PC2"]].mean().values
        vals, vecs = np.linalg.eigh(cov)
        order = vals.argsort()[::-1]
        vals, vecs = vals[order], vecs[:, order]
        theta = np.degrees(np.arctan2(*vecs[:, 0][::-1]))

        # Scale to 95% confidence region
        width, height = 2 * np.sqrt(chi2.ppf(0.95, 2)) * np.sqrt(vals)
        ellipse = Ellipse(
            xy=mean, width=width, height=height, angle=theta,
            edgecolor='black', facecolor='none',
            linestyle='--', linewidth=1.2
        )
        ax.add_patch(ellipse)

    # Final plot settings
    plt.title("PCA of Phenotypic Profiles (with 95% Ellipses)")
    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
    plt.legend(title="Clade")
    plt.grid(True)
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    ax.set_xlim(xlim[0] - 3, xlim[1] + 3)
    ax.set_ylim(ylim[0] - 3, ylim[1] + 3)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.savefig(output_file.replace(".pdf", ".png"), dpi=300)
    plt.close()


# -------------------- Section 5.1 – Clade-wise Phenotype Association Testing -------------------- #

def run_cladewise_stats(df, output_file):
    """
    Perform statistical test for each phenotype across clades.
    - Uses Fisher's exact test if only two clades; otherwise chi-squared test.
    - Applies Benjamini-Hochberg FDR correction.
    - Outputs CSV with p-values and adjusted q-values.
    """
    # Extract clade ID from strain names
    clades = [col.split("_")[0] for col in df.columns]

    # Ensure binary matrix (0/1 only)
    df_bin = df.fillna(0).astype(int)

    results = []

    for phenotype in df_bin.index:
        phenotype_data = df_bin.loc[phenotype]

        # Group phenotype presence/absence by clade
        counts_by_clade = {}
        for clade in set(clades):
            indices = [i for i, c in enumerate(clades) if c == clade]
            positives = phenotype_data.iloc[indices].sum()
            negatives = len(indices) - positives
            counts_by_clade[clade] = [positives, negatives]

        # Run appropriate test
        try:
            if len(counts_by_clade) == 2:
                # Fisher's exact test for 2x2 table
                _, p = fisher_exact(list(counts_by_clade.values()))
            else:
                # Chi-squared for >2 clades
                contingency = [counts_by_clade[c] for c in sorted(counts_by_clade)]
                _, p, _, _ = chi2_contingency(contingency)
        except ValueError:
            p = 1.0  # If test fails, treat as non-significant

        results.append((phenotype, p))

    # Multiple testing correction
    pvals = [r[1] for r in results]
    _, qvals, _, _ = multipletests(pvals, method='fdr_bh')

    # Save results
    pd.DataFrame({
        "Phenotype": [r[0] for r in results],
        "Raw_p": pvals,
        "FDR_q": qvals
    }).sort_values("FDR_q").to_csv(output_file, index=False)


# -------------------- Section 5.2 – Clade-wise Phenotype Association Testing -------------------- #

def plot_phenotype_correlation(df, prefix="phenotype_correlation"):
    """
    Compute Spearman correlation between phenotypes and visualize as a clustermap.
    Also saves:
      - Binary matrix used in the analysis
      - Raw correlation matrix as CSV
    """
    # Convert phenotype matrix to binary (0/1)
    df_bin = df.fillna(0).replace({"positive": 1, "negative": 0}).astype(int)
    df_bin.to_csv(f"{prefix}_binary_matrix.csv")
    print(f"[INFO] Binary matrix saved to '{prefix}_binary_matrix.csv'")

    # Remove phenotypes (rows) with no variation across strains
    df_var = df_bin.loc[df_bin.var(axis=1) > 0]
    if df_var.shape[0] < 2:
        print("[INFO] Too few variable phenotypes to compute correlation heatmap.")
        return

    # Transpose for phenotype-wise correlation, compute Spearman
    corr_df = df_var.T.corr(method="spearman")
    corr_df.to_csv(f"{prefix}_correlation_matrix.csv")
    print(f"[INFO] Correlation matrix saved to '{prefix}_correlation_matrix.csv'")

    # Remove rows/columns with NaN
    corr_df_clean = corr_df.dropna(axis=0).dropna(axis=1)
    if corr_df_clean.shape[0] < 2:
        print("[INFO] Correlation matrix has insufficient valid data for plotting.")
        return

    # Plot clustermap of correlations
    g = sns.clustermap(
        corr_df_clean,
        cmap="vlag",
        figsize=(20, 20),
        center=0,
        linewidths=0.1,
        xticklabels=True,
        yticklabels=True,
        cbar_kws={'label': 'Spearman correlation'}
    )
    g.fig.suptitle("Phenotype–Phenotype Correlation Heatmap", y=1.02)
    plt.tight_layout()

    # Save PDF and PNG
    g.savefig(f"{prefix}_clustermap.pdf")
    g.savefig(f"{prefix}_clustermap.png", dpi=300)
    plt.close()
    print(f"[INFO] Correlation clustermap saved as '{prefix}_clustermap.pdf'")


# -------------------- Section 5.3 – Random Forest Feature Importance Analysis -------------------- #

def run_random_forest_feature_importance(df, prefix="rf"):
    """
    Train a Random Forest classifier to predict clade based on phenotypic profiles.
    Outputs:
      - Feature importances ranked
      - Barplot of top 15 features
      - Classification report (precision, recall, F1-score)
      - Confusion matrix figure
    """
    # Transpose matrix: strains as rows, phenotypes as features
    df_bin = df.fillna(0).astype(int).T

    # Extract clade label (first part of strain name: A_Cneo_H99 → A)
    clade_labels = df_bin.index.str.split("_").str[0]
    X = df_bin.values
    y = clade_labels

    # Train Random Forest classifier
    clf = RandomForestClassifier(n_estimators=100, random_state=42)
    clf.fit(X, y)

    # Get and save feature importances
    importances = pd.Series(clf.feature_importances_, index=df_bin.columns).sort_values(ascending=False)
    importances.to_csv(f"{prefix}_rf_feature_importance.csv")
    print(f"[INFO] Random Forest feature importances saved to '{prefix}_rf_feature_importance.csv'")

    # Plot top 15 features
    plt.figure(figsize=(10, 6))
    importances.head(15).plot(kind='barh')
    plt.gca().invert_yaxis()
    plt.title("Top 15 Most Informative Phenotypes (Random Forest)")
    plt.xlabel("Feature Importance")
    plt.tight_layout()
    plt.savefig(f"{prefix}_rf_feature_importance.png", dpi=300)
    plt.savefig(f"{prefix}_rf_feature_importance.pdf")
    plt.close()

    # Classification metrics on training set
    y_pred = clf.predict(X)
    report = classification_report(y, y_pred)
    with open(f"{prefix}_rf_classification_report.txt", "w") as f:
        f.write(report)
    print(f"[INFO] Classification report saved to '{prefix}_rf_classification_report.txt'")

    # Confusion matrix
    cm = confusion_matrix(y, y_pred, labels=sorted(set(y)))
    disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=sorted(set(y)))
    fig, ax = plt.subplots(figsize=(6, 5))
    disp.plot(ax=ax, cmap="Blues")
    plt.title("Random Forest Confusion Matrix")
    plt.tight_layout()
    plt.savefig(f"{prefix}_rf_confusion_matrix.png", dpi=300)
    plt.savefig(f"{prefix}_rf_confusion_matrix.pdf")
    plt.close()



# -------------------- Section 5.4 – Barplot of Top Phenotypes by Clade -------------------- #

def plot_top_phenotype_clade_presence(df, feature_csv, prefix="rf", top_n=15):
    """
    Generate a barplot showing the average presence of the top N most informative phenotypes per clade.
    Uses Random Forest feature importance to select phenotypes.

    Inputs:
        - df: binary phenotype matrix (phenotypes × strains)
        - feature_csv: path to CSV with ranked features (from Random Forest)
        - top_n: number of top features to visualize
    """
    # Load top N phenotype features (most important)
    df_feat = pd.read_csv(feature_csv, index_col=0)
    importance_column = df_feat.columns[0]
    top_features = df_feat.sort_values(by=importance_column, ascending=False).head(top_n).index.tolist()

    # Subset binary matrix to top phenotypes
    df_bin = df.fillna(0).astype(int)
    df_top = df_bin.loc[top_features]

    # Decompose strain labels: MultiIndex (Clade, Strain)
    df_top.columns = pd.MultiIndex.from_tuples(
        [col.split("_", 1) for col in df_top.columns],
        names=["Clade", "Strain"]
    )

    # Average presence (0–1) of each phenotype per clade
    df_presence = df_top.T.groupby(level="Clade").mean().T  # shape: phenotype × clade

    # Assign distinct colors to phenotypes
    cmap = colormaps.get_cmap('tab20').resampled(top_n)
    colors = cmap.colors[:top_n]

    # Plot grouped barplot
    plt.figure(figsize=(12, 8))
    df_presence.plot(kind="bar", figsize=(12, 8), width=0.8, color=colors)

    plt.title(f"Presence Frequency of Top {top_n} Phenotypes by Clade")
    plt.ylabel("Mean Presence (0–1)")
    plt.xlabel("Clade")
    plt.xticks(rotation=90, ha='center', fontsize=9)
    plt.legend(title="Phenotype", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()

    plot_file = f"{prefix}_top{top_n}_phenotypes_by_clade_barplot.png"
    plt.savefig(plot_file, dpi=300)
    plt.savefig(plot_file.replace(".png", ".pdf"))
    plt.close()

    print(f"[INFO] Summary barplot saved to '{plot_file}' and PDF version.")


# -------------------- Section 5.5 – Duplicate Phenotype Detection -------------------- #

def find_duplicate_phenotypes(df, output_file="duplicate_phenotypes.csv"):
    """
    Identify phenotypes that have identical binary profiles across all strains.
    Useful to detect redundancy in the test set (e.g. synonymous or repeated substrates).

    Saves:
        - CSV file with columns: Representative, Duplicate
    """
    duplicates = {}  # map: representative → list of duplicates
    seen = {}        # map: binary profile → representative name

    for idx, row in df.iterrows():
        key = tuple(row.values)  # convert binary vector to a hashable tuple
        if key in seen:
            duplicates.setdefault(seen[key], []).append(idx)
        else:
            seen[key] = idx

    # Export duplicates if any were found
    if duplicates:
        with open(output_file, "w") as f:
            f.write("Representative,Duplicate\n")
            for rep, dups in duplicates.items():
                for dup in dups:
                    f.write(f"{rep},{dup}\n")
        print(f"[INFO] Found {sum(len(v) for v in duplicates.values())} duplicated phenotypes.")
        print(f"[INFO] Duplicate list saved to '{output_file}'")
    else:
        print("[INFO] No duplicate phenotypic profiles found.")



# -------------------- Main execution block -------------------- #

if __name__ == "__main__":
    import argparse

    # Set up command-line interface
    parser = argparse.ArgumentParser(description="Biolog phenotype analysis")
    parser.add_argument("input_file", help="Parsed phenotype matrix (.xlsx)")
    parser.add_argument("--with-plate-tags", action="store_true", help="Annotate tests with Biolog plate origin")
    parser.add_argument("--plate-file", help="TSV file for test-to-plate mapping (required with --with-plate-tags)")
    args = parser.parse_args()

    # Derive output prefix from input filename
    filepath = args.input_file
    prefix = os.path.splitext(os.path.basename(filepath))[0]

    # -------------------- Load & Preprocess -------------------- #
    df = load_and_clean_data(filepath)
    df.index = [greekify(name) for name in df.index]  # Convert names to Greek symbols (e.g., alpha → α)

    # -------------------- Optional Plate Annotation -------------------- #
    if args.with_plate_tags:
        if not args.plate_file:
            print("[ERROR] --plate-file is required with --with-plate-tags")
            sys.exit(1)
        df_to_plot = annotate_phenotype_with_plate_from_well_table(df, args.plate_file)
        annot_suffix = "_annot"
    else:
        df_to_plot = df
        annot_suffix = ""

    print(f"[INFO] Loaded {df.shape[0]} phenotypes across {df.shape[1]} strains")

    # -------------------- Visualization Outputs -------------------- #
    plot_clustermap_versions(df_to_plot, prefix=prefix + annot_suffix)
    plot_heatmap(df_to_plot, output_file=f"{prefix}{annot_suffix}_heatmap.png")
    plot_dendrogram(df_to_plot, output_file=f"{prefix}{annot_suffix}_dendrogram.pdf")
    plot_pca(df, output_file=f"{prefix}_pca.pdf")

    # -------------------- Statistical & ML Analysis -------------------- #
    run_cladewise_stats(df, output_file=f"{prefix}_cladewise_phenotype_stats.csv")
    plot_phenotype_correlation(df, prefix=f"{prefix}_phenotype")
    run_random_forest_feature_importance(df, prefix=prefix)
    plot_top_phenotype_clade_presence(df, feature_csv=f"{prefix}_rf_feature_importance.csv", prefix=prefix)

    # -------------------- Detect Redundant Phenotypes -------------------- #
    find_duplicate_phenotypes(df, output_file=f"{prefix}_duplicate_phenotypes.csv")


