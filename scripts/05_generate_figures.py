#!/usr/bin/env python3
"""This code will be to generate the plots for volcano and heatmap figures for the DGE analysis and this is for the poxvirus."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


REPO_ROOT = Path(__file__).resolve().parents[1]
PROCESSED_DIR = REPO_ROOT / "data" / "processed"
TABLE_DIR = REPO_ROOT / "results" / "tables"
FIGURE_DIR = REPO_ROOT / "results" / "figures"
KEY_LABELS = {
    "DHX29",
    "EIF1",
    "EIF1AX",
    "EIF2A",
    "EIF2S1",
    "EIF2S2",
    "EIF2S3",
    "EIF3A",
    "EIF3B",
    "EIF3C",
    "EIF3D",
    "EIF3E",
    "EIF3F",
    "EIF3G",
    "EIF3H",
    "EIF3I",
    "EIF3J",
    "EIF3K",
    "EIF3L",
    "EIF3M",
    "EIF4A1",
    "EIF4A2",
    "EIF4E",
    "EIF4G1",
    "EIF4G2",
    "EIF5",
    "EIF5A",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--padj", type=float, default=0.05, help="Adjusted p-value threshold for significance.")
    parser.add_argument("--top-n", type=int, default=30, help="Number of translation factors for the heatmap.")
    return parser.parse_args()


def classify_row(row: pd.Series, padj_threshold: float) -> str:
    if pd.isna(row["padj"]) or row["padj"] >= padj_threshold:
        return "not significant"
    if row["log2FoldChange"] > 0:
        return "upregulated"
    if row["log2FoldChange"] < 0:
        return "downregulated"
    return "not significant"


def safe_neg_log10(values: pd.Series) -> pd.Series:
    finite_positive = values[(values > 0) & np.isfinite(values)]
    floor = finite_positive.min() / 10 if not finite_positive.empty else 1e-300
    return -np.log10(values.fillna(1.0).clip(lower=floor))


def generate_volcano(results: pd.DataFrame, padj: float) -> list[str]:
    plot_df = results.copy()
    plot_df["neg_log10_padj"] = safe_neg_log10(plot_df["padj"])
    plot_df["class"] = plot_df.apply(classify_row, axis=1, padj_threshold=padj)

    palette = {
        "not significant": "#808080",
        "upregulated": "#C43C39",
        "downregulated": "#2F6DB3",
    }

    sns.set_theme(style="whitegrid", context="paper")
    plt.figure(figsize=(9, 6.5))
    ax = sns.scatterplot(
        data=plot_df,
        x="log2FoldChange",
        y="neg_log10_padj",
        hue="class",
        palette=palette,
        s=14,
        linewidth=0,
        alpha=0.72,
    )
    ax.axhline(-math.log10(padj), color="#2B2B2B", linestyle="--", linewidth=0.9)
    ax.axvline(0, color="#2B2B2B", linestyle="-", linewidth=0.7)
    ax.set_xlabel("log2 fold change (VacV vs mock)")
    ax.set_ylabel("-log10 adjusted p-value")
    ax.set_title("GSE278320 Host Differential Expression")

    labeled: list[str] = []
    label_df = plot_df[plot_df["gene_symbol"].isin(KEY_LABELS)].copy()
    label_df = label_df.sort_values(["padj", "gene_symbol"], na_position="last")
    for _, row in label_df.iterrows():
        if pd.isna(row["log2FoldChange"]) or pd.isna(row["neg_log10_padj"]):
            continue
        ax.text(
            row["log2FoldChange"],
            row["neg_log10_padj"],
            str(row["gene_symbol"]),
            fontsize=7,
            ha="left",
            va="bottom",
        )
        labeled.append(str(row["gene_symbol"]))

    ax.legend(title="", frameon=True, loc="upper right")
    plt.tight_layout()
    fig = ax.get_figure()
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(FIGURE_DIR / f"volcano_plot.{suffix}", dpi=300, bbox_inches="tight")
    plt.close()
    return labeled


def generate_heatmap(translation: pd.DataFrame, normalized_counts: pd.DataFrame, metadata: pd.DataFrame, top_n: int) -> list[str]:
    candidates = translation[translation["padj"].notna()].sort_values(["padj", "log2FoldChange"]).head(top_n)
    if candidates.empty:
        candidates = translation.sort_values("baseMean", ascending=False).head(top_n)
    if candidates.empty:
        raise SystemExit("No translation-factor rows are available for heatmap generation.")

    gene_ids = [gene_id for gene_id in candidates["gene_id"].tolist() if gene_id in normalized_counts.columns]
    if not gene_ids:
        raise SystemExit("No selected translation-factor genes were found in normalized_counts.csv.")

    labels = candidates.set_index("gene_id").loc[gene_ids, "gene_symbol"].fillna(pd.Series(gene_ids, index=gene_ids))
    heatmap_data = np.log2(normalized_counts[gene_ids] + 1).T
    heatmap_data.index = [f"{labels.loc[gene_id]} ({gene_id})" for gene_id in gene_ids]

    sample_labels = metadata["infection"].astype(str) + "_rep" + metadata["replicate"].astype(str)
    heatmap_data.columns = [sample_labels.loc[sample_id] for sample_id in heatmap_data.columns]

    sns.set_theme(style="white", context="paper")
    cluster_grid = sns.clustermap(
        heatmap_data,
        cmap="vlag",
        z_score=0,
        figsize=(10, 9),
        col_cluster=False,
        yticklabels=True,
        xticklabels=True,
        cbar_kws={"label": "row z-score of log2 normalized counts"},
    )
    cluster_grid.fig.suptitle("Translation-Factor Expression Across GSE278320 Total-Fraction Samples", y=1.02)
    cluster_grid.ax_heatmap.set_xlabel("")
    cluster_grid.ax_heatmap.set_ylabel("")
    for suffix in ("png", "pdf", "svg"):
        cluster_grid.fig.savefig(FIGURE_DIR / f"translation_heatmap.{suffix}", dpi=300, bbox_inches="tight")
    plt.close(cluster_grid.fig)
    return labels.tolist()


def main() -> None:
    args = parse_args()
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    results = pd.read_csv(TABLE_DIR / "dge_results_full.csv", dtype={"gene_id": "string", "gene_symbol": "string"})
    translation = pd.read_csv(TABLE_DIR / "translation_factors_all.csv", dtype={"gene_id": "string", "gene_symbol": "string"})
    normalized_counts = pd.read_csv(PROCESSED_DIR / "normalized_counts.csv", index_col="sample_id")
    metadata = pd.read_csv(PROCESSED_DIR / "model_ready_metadata.csv", index_col="sample_id")

    labeled = generate_volcano(results, padj=args.padj)
    heatmap_genes = generate_heatmap(translation, normalized_counts, metadata, top_n=args.top_n)

    available_targets = sorted(set(results["gene_symbol"].dropna()) & KEY_LABELS)
    unlabeled = sorted(set(available_targets) - set(labeled))
    print(f"Volcano labels: {', '.join(labeled) if labeled else 'none'}")
    if unlabeled:
        print(f"Key genes present but not labeled due to missing plot coordinates: {', '.join(unlabeled)}")
    print(f"Heatmap genes: {', '.join(heatmap_genes[:10])}{' ...' if len(heatmap_genes) > 10 else ''}")
    print(f"Wrote {FIGURE_DIR / 'volcano_plot.png'} and {FIGURE_DIR / 'translation_heatmap.png'}.")


if __name__ == "__main__":
    main()
