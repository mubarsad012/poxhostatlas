#!/usr/bin/env python3
"""This is for creating more of an advanced synthesis figures for expanded cross-dataset poxvirus analysis.

What the below code will be doing is generating publication panels not covered by the core pipeline:
  -Host genetic-background robustness heatmap (GSE288000 five-way panel)
  - Vaccinia temporal dynamics for conserved translation factors (GSE284044)
  -Pan-meta summary for directional concordance for helicase/eIF modules
"""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

REPO_ROOT = Path(__file__).resolve().parents[1]
EXPANDED_DIR = REPO_ROOT / "results" / "expanded"
META_DIR = REPO_ROOT / "results" / "meta_analysis"
OUT_DIR = REPO_ROOT / "results" / "synthesis"
FIG_DIR = OUT_DIR / "figures"

TARGET_PATTERN = re.compile(r"^(DHX|DDX|EIF|RPS|RPL)", re.IGNORECASE)
GSE288000_BACKGROUNDS = ["NTC", "N4BP1", "ZC3H12A", "TRIM25", "ZC3H12A_N4BP1"]
TIMEPOINTS = [2, 6, 24]
FOCAL_GENES = ["DHX29", "DHX15", "DDX21", "EIF4B", "EIF3L", "RPL10", "RPS4X"]


def host_background_heatmap() -> None:
    path = EXPANDED_DIR / "GSE288000_host_background_robustness.csv"
    if not path.exists():
        return
    frame = pd.read_csv(path, index_col=0)
    frame.index = frame.index.astype(str).str.upper()
    tf = frame[frame.index.str.match(TARGET_PATTERN)]
    tf = tf.dropna(subset=["direction_consistency"], how="all")
    tf = tf[tf["direction_consistency"] >= 0.8].sort_values("mean_log2FoldChange", ascending=False).head(35)
    plot = tf[GSE288000_BACKGROUNDS].astype(float)
    sns.set_theme(style="white", context="paper")
    fig, ax = plt.subplots(figsize=(7.5, max(4.5, len(plot) * 0.28)))
    sns.heatmap(
        plot,
        cmap="RdBu_r",
        center=0,
        vmin=-1.5,
        vmax=1.5,
        linewidths=0.3,
        linecolor="#E8E8E8",
        ax=ax,
        cbar_kws={"label": "log2 fold change (M003 vs mCherry)"},
    )
    ax.set_title("Translation-factor remodeling across GSE288000 host backgrounds")
    ax.set_xlabel("Host genetic background")
    ax.set_ylabel("")
    plt.setp(ax.get_xticklabels(), rotation=30, ha="right")
    fig.tight_layout()
    for suf in ("png", "pdf", "svg"):
        fig.savefig(FIG_DIR / f"host_background_robustness_heatmap.{suf}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def temporal_dynamics_panel() -> None:
    rows = []
    for hpi in TIMEPOINTS:
        path = EXPANDED_DIR / f"GSE284044_{hpi}hpi_dge.csv"
        if not path.exists():
            continue
        d = pd.read_csv(path)
        d["gene_symbol"] = d["gene_symbol"].astype(str).str.upper()
        d["hpi"] = hpi
        rows.append(d[d["gene_symbol"].isin(FOCAL_GENES)][["gene_symbol", "hpi", "log2FoldChange", "padj"]])
    if not rows:
        return
    long = pd.concat(rows, ignore_index=True)
    sns.set_theme(style="whitegrid", context="paper")
    fig, ax = plt.subplots(figsize=(8, 4.5))
    palette = sns.color_palette("tab10", n_colors=len(FOCAL_GENES))
    for (gene, sub), color in zip(long.groupby("gene_symbol"), palette):
        sub = sub.sort_values("hpi")
        ax.plot(sub["hpi"], sub["log2FoldChange"], marker="o", lw=2, label=gene, color=color)
        sig = sub["padj"] < 0.05
        ax.scatter(sub.loc[sig, "hpi"], sub.loc[sig, "log2FoldChange"], s=120, facecolors="none",
                   edgecolors=color, linewidths=2, zorder=5)
    ax.axhline(0, color="#333", lw=0.8, ls="--")
    ax.set_xticks(TIMEPOINTS)
    ax.set_xlabel("Hours post Vaccinia infection (GSE284044, Vero)")
    ax.set_ylabel("log2 fold change vs mock")
    ax.set_title("Temporal translation-factor dynamics (open circles = FDR < 0.05)")
    ax.legend(frameon=False, fontsize=8, ncol=2, loc="best")
    fig.tight_layout()
    for suf in ("png", "pdf", "svg"):
        fig.savefig(FIG_DIR / f"vaccinia_temporal_dynamics.{suf}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def module_concordance_bar() -> None:
    pan_path = META_DIR / "meta_pan_poxvirus_translation_factors.csv"
    if not pan_path.exists():
        return
    meta = pd.read_csv(pan_path, index_col=0)
    meta.index = meta.index.astype(str).str.upper()
    modules = {
        "DHX helicases": meta.index.str.match(r"^DHX"),
        "DDX helicases": meta.index.str.match(r"^DDX"),
        "EIF factors": meta.index.str.match(r"^EIF"),
        "RPL ribosomal": meta.index.str.match(r"^RPL"),
        "RPS ribosomal": meta.index.str.match(r"^RPS"),
    }
    records = []
    for name, mask in modules.items():
        sub = meta[mask & meta["k_studies"].ge(3)]
        if sub.empty:
            continue
        records.append({
            "module": name,
            "n_genes": len(sub),
            "pct_up_meta": 100 * (sub["pooled_log2FoldChange"] > 0).mean(),
            "pct_concordant": 100 * sub["direction_concordance"].mean(),
            "median_pooled_lfc": sub["pooled_log2FoldChange"].median(),
            "sig_meta_FDR05": int((sub["meta_FDR"] < 0.05).sum()),
        })
    summary = pd.DataFrame(records)
    summary.to_csv(OUT_DIR / "module_concordance_summary.csv", index=False)
    sns.set_theme(style="whitegrid", context="paper")
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.8))
    x = np.arange(len(summary))
    axes[0].bar(x, summary["pct_concordant"], color="#4E79A7", edgecolor="white")
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(summary["module"], rotation=25, ha="right")
    axes[0].set_ylabel("% directional concordance\n(across 4 studies)")
    axes[0].set_ylim(0, 105)
    axes[0].set_title("Cross-study direction consistency")
    colors = ["#B94E48" if v > 0 else "#3E6FB6" for v in summary["median_pooled_lfc"]]
    axes[1].bar(x, summary["median_pooled_lfc"], color=colors, edgecolor="white")
    axes[1].axhline(0, color="#333", lw=0.8)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(summary["module"], rotation=25, ha="right")
    axes[1].set_ylabel("Median pooled log2FC")
    axes[1].set_title("Pan-poxvirus random-effects effect")
    fig.suptitle("Translation-module meta-analysis synthesis", y=1.02, fontsize=11)
    fig.tight_layout()
    for suf in ("png", "pdf", "svg"):
        fig.savefig(FIG_DIR / f"module_concordance_summary.{suf}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    host_background_heatmap()
    temporal_dynamics_panel()
    module_concordance_bar()
    print(f"Advanced synthesis figures written to {FIG_DIR}")


if __name__ == "__main__":
    main()
