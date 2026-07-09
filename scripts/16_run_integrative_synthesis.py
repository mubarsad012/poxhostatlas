#!/usr/bin/env python3
"""Integrative evidence synthesis: cross-study concordance, external validation, pipeline figure.

Extends the meta-analysis with:
  - GSE185520 independent Vaccinia host-shutoff directional validation (external XLSX)
  - Cross-study Spearman rank correlation of translation-factor effect sizes
  - Composite evidence score ranking novel multi-context hits
  - Extended evidence heatmap and pipeline schematic for manuscript embedding
"""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

REPO_ROOT = Path(__file__).resolve().parents[1]
TABLE_DIR = REPO_ROOT / "results" / "tables"
META_DIR = REPO_ROOT / "results" / "meta_analysis"
EXPANDED_DIR = REPO_ROOT / "results" / "expanded"
EXTERNAL_RAW = REPO_ROOT / "data" / "external" / "expansion_raw"
OUT_DIR = REPO_ROOT / "results" / "synthesis"
FIG_DIR = OUT_DIR / "figures"

TARGET_PATTERN = re.compile(r"^(DHX|DDX|EIF|RPS|RPL)", re.IGNORECASE)

GSE185520_CONTRASTS = {
    "GSE185520 WT 6hr": "GSE185520_MockvsWT_6hr.xlsx",
    "GSE185520 WT 18hr": "GSE185520_MockvsWT_18hr.xlsx",
    "GSE185520 D9muD10mu 18hr": "GSE185520_MockvsD9muD10mu_18hr.xlsx",
}

CORE_CONTEXTS = {
    "GSE278320 VacV": TABLE_DIR / "dge_results_full.csv",
    "GSE284044 6hpi": EXPANDED_DIR / "GSE284044_6hpi_dge.csv",
    "GSE287860 M003": REPO_ROOT / "results" / "meta" / "GSE287860_dge_results.csv",
    "GSE288000 NTC": REPO_ROOT / "results" / "meta" / "GSE288000_NTC_dge_results.csv",
}


def load_harmonized(path: Path, gene_col: str = "gene_symbol") -> pd.DataFrame:
    frame = pd.read_csv(path)
    frame[gene_col] = frame[gene_col].astype(str).str.upper()
    return frame[[gene_col, "log2FoldChange"]].rename(
        columns={gene_col: "gene_symbol", "log2FoldChange": "log2fc"}
    )


def load_gse185520(filename: str) -> pd.DataFrame:
    path = EXTERNAL_RAW / filename
    if not path.exists():
        return pd.DataFrame(columns=["gene_symbol", "log2fc"])
    frame = pd.read_excel(path)
    frame["gene_symbol"] = frame["gene"].astype(str).str.upper()
    frame = frame.rename(columns={"log2(fold_change)": "log2fc"})
    return frame[["gene_symbol", "log2fc"]].dropna(subset=["log2fc"])


def build_extended_context_table() -> pd.DataFrame:
    rows = []
    for ctx, path in CORE_CONTEXTS.items():
        if not path.exists():
            continue
        sub = load_harmonized(path)
        sub["context"] = ctx
        rows.append(sub)
    for ctx, fname in GSE185520_CONTRASTS.items():
        sub = load_gse185520(fname)
        if sub.empty:
            continue
        sub["context"] = ctx
        rows.append(sub)
    return pd.concat(rows, ignore_index=True)


def rank_correlation_matrix(long: pd.DataFrame) -> pd.DataFrame:
    tf = long[long["gene_symbol"].str.match(TARGET_PATTERN)]
    pivot = tf.pivot_table(index="gene_symbol", columns="context", values="log2fc", aggfunc="mean")
    corr = pivot.corr(method="spearman")
    corr.to_csv(OUT_DIR / "cross_study_spearman_correlation.csv")
    return corr


def plot_rank_correlation(corr: pd.DataFrame) -> None:
    if corr.empty:
        return
    sns.set_theme(style="white", context="paper")
    fig, ax = plt.subplots(figsize=(8.5, 6.5))
    sns.heatmap(
        corr,
        annot=True,
        fmt=".2f",
        cmap="RdYlBu_r",
        center=0,
        vmin=-1,
        vmax=1,
        square=True,
        linewidths=0.5,
        ax=ax,
        cbar_kws={"label": "Spearman rho (translation-factor log2FC ranks)"},
    )
    ax.set_title("Cross-study rank concordance of translation-factor remodeling")
    plt.setp(ax.get_xticklabels(), rotation=35, ha="right", fontsize=8)
    plt.setp(ax.get_yticklabels(), fontsize=8)
    fig.tight_layout()
    for suf in ("png", "pdf", "svg"):
        fig.savefig(FIG_DIR / f"cross_study_rank_correlation.{suf}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def composite_evidence_score(long: pd.DataFrame) -> pd.DataFrame:
    tf = long[long["gene_symbol"].str.match(TARGET_PATTERN)].copy()
    pivot = tf.pivot_table(index="gene_symbol", columns="context", values="log2fc", aggfunc="mean")
    meta = pd.read_csv(META_DIR / "meta_pan_poxvirus_translation_factors.csv", index_col=0)
    meta.index = meta.index.astype(str).str.upper()

    records = []
    for gene in pivot.index:
        row = pivot.loc[gene].dropna()
        if row.empty:
            continue
        n_ctx = len(row)
        n_up = int((row > 0).sum())
        n_down = int((row < 0).sum())
        direction_consistency = max(n_up, n_down) / n_ctx
        meta_row = meta.loc[gene] if gene in meta.index else None
        records.append(
            {
                "gene_symbol": gene,
                "n_contexts": n_ctx,
                "mean_log2fc": row.mean(),
                "direction_consistency": direction_consistency,
                "n_positive_contexts": n_up,
                "n_negative_contexts": n_down,
                "pooled_log2fc_meta": float(meta_row["pooled_log2FoldChange"]) if meta_row is not None else np.nan,
                "meta_FDR": float(meta_row["meta_FDR"]) if meta_row is not None else np.nan,
                "meta_k_studies": int(meta_row["k_studies"]) if meta_row is not None else 0,
                **{f"log2fc_{c.replace(' ', '_')}": row.get(c, np.nan) for c in pivot.columns},
            }
        )
    score = pd.DataFrame(records)
    score["meta_sig"] = score["meta_FDR"].lt(0.05)
    score["evidence_score"] = (
        score["direction_consistency"] * 3
        + score["n_contexts"] * 0.5
        + score["meta_sig"].astype(float) * 4
        + np.clip(score["pooled_log2fc_meta"].abs(), 0, 2) * 0.5
    )
    score = score.sort_values(["evidence_score", "meta_FDR"], ascending=[False, True])
    score.to_csv(OUT_DIR / "composite_evidence_score.csv", index=False)
    return score


def plot_extended_evidence_heatmap(long: pd.DataFrame, top_n: int = 35) -> None:
    score = pd.read_csv(OUT_DIR / "composite_evidence_score.csv")
    genes = score.head(top_n)["gene_symbol"].tolist()
    pivot = long[long["gene_symbol"].isin(genes)].pivot_table(
        index="gene_symbol", columns="context", values="log2fc", aggfunc="mean"
    )
    pivot = pivot.reindex(genes)
    col_order = list(CORE_CONTEXTS.keys()) + list(GSE185520_CONTRASTS.keys())
    pivot = pivot.reindex(columns=[c for c in col_order if c in pivot.columns])
    if pivot.empty:
        return
    sns.set_theme(style="white", context="paper")
    fig, ax = plt.subplots(figsize=(11, max(5, len(pivot) * 0.32)))
    sns.heatmap(
        pivot,
        cmap="vlag",
        center=0,
        vmin=-2.5,
        vmax=2.5,
        linewidths=0.3,
        linecolor="#E6E6E6",
        ax=ax,
        cbar_kws={"label": "log2 fold change"},
    )
    ax.set_title("Extended translation-factor evidence across core + GSE185520 contexts")
    ax.set_xlabel("")
    ax.set_ylabel("")
    plt.setp(ax.get_xticklabels(), rotation=40, ha="right", fontsize=7.5)
    fig.tight_layout()
    for suf in ("png", "pdf", "svg"):
        fig.savefig(FIG_DIR / f"extended_evidence_heatmap.{suf}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_pipeline_overview() -> None:
    labels = [
        "GEO ingestion\n(4 core + GSE185520)",
        "PyDESeq2 DGE\n+ QC manifest",
        "GENCODE host\ncuration",
        "10 harmonized\ncontrasts",
        "Random-effects\nmeta-analysis",
        "Robustness +\ntemporal panels",
        "STRING network\n+ evidence score",
        "Manuscript +\nrelease archive",
    ]
    colors = ["#E8EEF7", "#EAF4EC", "#FFF2CC", "#FCE4D6", "#EADCF8", "#DDEAF6", "#F4E4D6", "#E2F0F3"]
    fig, ax = plt.subplots(figsize=(13.5, 2.8))
    ax.set_axis_off()
    xs = np.linspace(0.06, 0.94, len(labels))
    for x, label, color in zip(xs, labels, colors, strict=True):
        ax.add_patch(plt.Rectangle((x - 0.055, 0.35), 0.11, 0.35, fc=color, ec="#555", lw=0.8))
        ax.text(x, 0.525, label, ha="center", va="center", fontsize=8.5, wrap=True)
    for i in range(len(xs) - 1):
        ax.annotate("", xy=(xs[i + 1] - 0.06, 0.525), xytext=(xs[i] + 0.06, 0.525),
                    arrowprops=dict(arrowstyle="->", color="#444", lw=1.2))
    ax.set_xlim(0, 1)
    ax.set_ylim(0.2, 0.85)
    ax.set_title("Cross-poxvirus translation-factor analysis pipeline", fontsize=11, pad=8)
    fig.tight_layout()
    for suf in ("png", "pdf", "svg"):
        fig.savefig(FIG_DIR / f"pipeline_overview.{suf}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def novel_findings_summary(score: pd.DataFrame) -> pd.DataFrame:
    dhx = score[score["gene_symbol"].str.match(r"^DHX")]
    helicase_up = score[
        score["gene_symbol"].str.match(r"^(DHX|DDX)")
        & score["meta_FDR"].lt(0.05)
        & score["pooled_log2fc_meta"].gt(0)
    ]
    eif_down = score[
        score["gene_symbol"].str.match(r"^EIF")
        & score["meta_FDR"].lt(0.05)
        & score["pooled_log2fc_meta"].lt(0)
    ]
    rows = [
        {"finding": "Pan-poxvirus meta-significant translation factors", "count": int(score["meta_sig"].sum())},
        {"finding": "Meta-significant DHX/DDX helicases (up)", "count": len(helicase_up)},
        {"finding": "Meta-significant EIF factors (down)", "count": len(eif_down)},
        {"finding": "DHX genes with >=80% direction consistency (all contexts)", "count": int((dhx["direction_consistency"] >= 0.8).sum())},
        {"finding": "Top composite evidence score", "count": float(score["evidence_score"].iloc[0]) if len(score) else np.nan},
        {"finding": "Top composite evidence gene", "count": score["gene_symbol"].iloc[0] if len(score) else ""},
    ]
    summary = pd.DataFrame(rows)
    summary.to_csv(OUT_DIR / "novel_findings_summary.csv", index=False)
    return summary


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    long = build_extended_context_table()
    long.to_csv(OUT_DIR / "extended_context_effects_long.csv", index=False)

    corr = rank_correlation_matrix(long)
    plot_rank_correlation(corr)

    score = composite_evidence_score(long)
    plot_extended_evidence_heatmap(long)
    plot_pipeline_overview()
    summary = novel_findings_summary(score)

    print("Integrative synthesis complete.")
    print(summary.to_string(index=False))
    dhx29 = score[score["gene_symbol"] == "DHX29"]
    if not dhx29.empty:
        r = dhx29.iloc[0]
        print(f"\nDHX29: evidence_score={r['evidence_score']:.2f} "
              f"contexts={int(r['n_contexts'])} consistency={r['direction_consistency']:.0%}")


if __name__ == "__main__":
    main()
