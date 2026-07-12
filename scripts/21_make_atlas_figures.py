#!/usr/bin/env python3
"""This is for rendering the six polished PoxHostAtlas main figures as composite PNG/PDF panels."""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import FancyArrowPatch, Rectangle

REPO_ROOT = Path(__file__).resolve().parents[1]
META_DIR = REPO_ROOT / "results" / "meta_analysis"
ML_DIR = REPO_ROOT / "results" / "ml"
NET_DIR = REPO_ROOT / "results" / "network"
SYN_DIR = REPO_ROOT / "results" / "synthesis"
MECH_DIR = REPO_ROOT / "results" / "mechanistic"
DOCS = REPO_ROOT / "docs"
FIG_DIR = REPO_ROOT / "results" / "figures" / "atlas"

TARGET_PATTERN = re.compile(r"^(DHX|DDX|EIF|RPS|RPL)", re.IGNORECASE)
MODULE_COLORS = {
    "DHX/DDX helicase": "#C0392B", "eIF3 complex": "#E67E22", "eIF4/initiation": "#F1C40F",
    "other eIF": "#16A085", "40S ribosomal": "#2980B9", "60S ribosomal": "#8E44AD", "other": "#95A5A6",
}


def save(fig, stub):
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    for suf in ("png", "pdf"):
        fig.savefig(FIG_DIR / f"{stub}.{suf}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def fig1_atlas():
    reg = pd.read_csv(DOCS / "dataset_registry.csv")
    fig = plt.figure(figsize=(14, 7.2))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.05, 0.95], hspace=0.38, wspace=0.22)

    # A. PRISMA-style flow (tighter)
    axA = fig.add_subplot(gs[0, 0]); axA.set_axis_off(); axA.set_title("A  Systematic dataset discovery", loc="left", fontweight="bold")
    nA = int((reg["tier"] == "A").sum()); nB = int((reg["tier"] == "B").sum()); nC = int((reg["tier"] == "C").sum())
    steps = [
        (f"GEO / SRA / PubMed search\n(poxvirus host transcriptomes)", "#D6EAF8"),
        (f"Candidates screened (n={len(reg)})", "#D5F5E3"),
        (f"Tier A reanalyzed (n={nA})  +  Tier B validation (n={nB})", "#FCF3CF"),
        (f"Tier C contextual / excluded (n={nC})", "#FADBD8"),
    ]
    y = 0.95
    for txt, col in steps:
        axA.add_patch(Rectangle((0.04, y - 0.17), 0.92, 0.155, fc=col, ec="#555"))
        axA.text(0.5, y - 0.092, txt, ha="center", va="center", fontsize=9)
        if y < 0.95:
            axA.annotate("", xy=(0.5, y + 0.0), xytext=(0.5, y + 0.05),
                         arrowprops=dict(arrowstyle="-|>", color="#444"))
        y -= 0.245
    axA.set_xlim(0, 1); axA.set_ylim(0.0, 1)

    # B. Samples per Tier-A study
    axB = fig.add_subplot(gs[0, 1])
    a = reg[reg["tier"] == "A"]
    axB.barh(a["accession"], a["n_samples"], color="#4E79A7", ec="white")
    for i, (acc, n) in enumerate(zip(a["accession"], a["n_samples"])):
        axB.text(n + 0.5, i, str(int(n)), va="center", fontsize=8)
    axB.set_xlabel("samples in series"); axB.set_title("B  Tier-A study scale", loc="left", fontweight="bold")
    axB.invert_yaxis()

    # C. Transparent tiering table
    axC = fig.add_subplot(gs[1, 0]); axC.set_axis_off()
    axC.set_title("C  Data tiering and usage", loc="left", fontweight="bold")
    tier_rows = [
        ["Study", "Tier", "Raw\ncounts", "Re-run\nDGE", "Meta", "ML"],
        ["GSE278320", "A", "yes", "yes", "yes", "yes"],
        ["GSE287860", "A", "yes", "yes", "yes", "yes"],
        ["GSE288000", "A", "yes", "yes", "yes", "yes"],
        ["GSE284044", "A", "DE tbl", "no*", "yes", "no"],
        ["GSE185520", "B", "no", "no", "no", "valid."],
    ]
    tbl = axC.table(cellText=tier_rows[1:], colLabels=tier_rows[0], loc="center", cellLoc="center")
    tbl.auto_set_font_size(False); tbl.set_fontsize(8); tbl.scale(1, 1.5)
    for j in range(len(tier_rows[0])):
        tbl[0, j].set_facecolor("#34495E"); tbl[0, j].get_text().set_color("white")
    axC.text(0.0, -0.02, "*author DESeq2 time-course tables harmonized; not re-run from raw counts",
             transform=axC.transAxes, fontsize=7, style="italic")

    # D. Pipeline strip
    axD = fig.add_subplot(gs[1, 1]); axD.set_axis_off(); axD.set_title("D  Analysis pipeline", loc="left", fontweight="bold")
    labels = ["DGE\n(PyDESeq2)", "Meta + LOSO\n+ het. classes", "LODO ML\n+ ablation + null", "Network +\nevidence score"]
    cols = ["#D5F5E3", "#FCF3CF", "#EBDEF0", "#D6EAF8"]
    ys = [0.8, 0.57, 0.34, 0.11]
    for yy, lab, c in zip(ys, labels, cols):
        axD.add_patch(Rectangle((0.1, yy - 0.08), 0.8, 0.15, fc=c, ec="#555"))
        axD.text(0.5, yy - 0.005, lab, ha="center", va="center", fontsize=8.5)
        if yy > 0.11:
            axD.annotate("", xy=(0.5, yy - 0.09), xytext=(0.5, yy - 0.14),
                         arrowprops=dict(arrowstyle="-|>", color="#444"))
    axD.set_xlim(0, 1); axD.set_ylim(0, 0.95)
    fig.suptitle("Figure 1. PoxHostAtlas: systematic public poxvirus host-response atlas, tiering, and pipeline", fontsize=13, fontweight="bold")
    save(fig, "Figure1_atlas")


def _forest(ax, gene, ctx):
    sub = ctx[ctx["gene_symbol"] == gene].dropna(subset=["log2FoldChange", "lfcSE"]).copy()
    sub = sub.iloc[::-1]
    ys = np.arange(len(sub))
    lo = sub["log2FoldChange"] - 1.96 * sub["lfcSE"]
    hi = sub["log2FoldChange"] + 1.96 * sub["lfcSE"]
    cols = ["#C0392B" if ("VacV" in c or "hpi" in c) else "#4E79A7" for c in sub["context"]]
    for y, (_, r), c, l, h in zip(ys, sub.iterrows(), cols, lo, hi):
        ax.plot([l, h], [y, y], color=c, lw=1.6)
        sig = pd.notna(r["padj"]) and r["padj"] < 0.05
        ax.scatter([r["log2FoldChange"]], [y], s=60 if sig else 38, color=c,
                   ec="black" if sig else "none", lw=0.8, zorder=3)
    ax.axvline(0, ls="--", color="#333", lw=0.8)
    ax.set_yticks(ys); ax.set_yticklabels(sub["context"], fontsize=7)
    ax.set_xlabel("log2FC (95% CI)"); ax.set_title(gene, fontweight="bold", fontsize=10)


def fig2_meta():
    pan = pd.read_csv(META_DIR / "meta_pan_poxvirus_full.csv")
    pan["gene"] = pan["gene_symbol"].astype(str).str.upper()
    ctx = pd.read_csv(META_DIR / "all_context_effects_long.csv")
    ctx["gene_symbol"] = ctx["gene_symbol"].astype(str).str.upper()
    cls = pd.read_csv(META_DIR / "heterogeneity_classification.csv")

    fig = plt.figure(figsize=(15, 10))
    gs = fig.add_gridspec(2, 3, hspace=0.30, wspace=0.30)

    # A volcano
    axA = fig.add_subplot(gs[0, :2])
    p = pan[pan["k_studies"] >= 3].dropna(subset=["meta_FDR", "pooled_log2FoldChange"]).copy()
    p["ml10"] = -np.log10(p["meta_FDR"].clip(lower=1e-300))
    tf = p[p["is_translation_factor"] == True]
    axA.scatter(p["pooled_log2FoldChange"], p["ml10"], s=8, c="#C9CDD2", alpha=0.4, lw=0)
    sig_tf = tf[tf["meta_FDR"] < 0.05]
    axA.scatter(sig_tf["pooled_log2FoldChange"], sig_tf["ml10"], s=26, c="#6B4EA0", alpha=0.85, lw=0, label="sig. translation factor")
    axA.axhline(-np.log10(0.05), ls="--", color="#333", lw=0.8)
    axA.axvline(0, color="#333", lw=0.8)
    for _, r in sig_tf.sort_values("meta_FDR").head(12).iterrows():
        axA.text(r["pooled_log2FoldChange"], r["ml10"], r["gene"], fontsize=7)
    axA.set_xlabel("pooled log2FC (random effects)"); axA.set_ylabel("-log10 meta-FDR")
    axA.set_title("A  Pan-poxvirus meta-analysis volcano", loc="left", fontweight="bold")
    axA.legend(frameon=False, fontsize=8)

    # B het class scatter
    axB = fig.add_subplot(gs[0, 2])
    palette = {"I_conserved_low_het": "#2ECC71", "II_conserved_high_het": "#F39C12",
               "III_context_dependent": "#3498DB", "IV_non_reproducible": "#E74C3C"}
    for cl, grp in cls.groupby("class"):
        axB.scatter(grp["direction_concordance"], grp["I2"], s=22, c=palette.get(cl, "#999"),
                    label=cl.split("_")[0], alpha=0.8, lw=0)
    for g in ["DHX15", "DHX29"]:
        r = cls[cls["gene"] == g]
        if not r.empty:
            axB.annotate(g, (r["direction_concordance"].iloc[0], r["I2"].iloc[0]),
                         fontsize=8, fontweight="bold")
    axB.set_xlabel("direction concordance"); axB.set_ylabel("I^2 (%)")
    axB.set_title("B  Heterogeneity classes", loc="left", fontweight="bold")
    axB.legend(frameon=False, fontsize=7, title_fontsize=7)

    # C/D forests
    _forest(fig.add_subplot(gs[1, 0]), "DHX15", ctx)
    _forest(fig.add_subplot(gs[1, 1]), "DHX29", ctx)

    # E module concordance
    axE = fig.add_subplot(gs[1, 2])
    mod = pd.read_csv(SYN_DIR / "module_concordance_summary.csv")
    colors = ["#C0392B" if v > 0 else "#2980B9" for v in mod["median_pooled_lfc"]]
    axE.barh(mod["module"], mod["median_pooled_lfc"], color=colors, ec="white")
    axE.axvline(0, color="#333", lw=0.8)
    axE.set_xlabel("median pooled log2FC"); axE.set_title("C  Module remodeling", loc="left", fontweight="bold")
    fig.suptitle("Figure 2. Cross-study meta-analysis and heterogeneity classification", fontsize=13, fontweight="bold")
    save(fig, "Figure2_meta")


def fig3_modules():
    fig = plt.figure(figsize=(13, 5.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.2, 1], wspace=0.3)
    # rank concordance
    axA = fig.add_subplot(gs[0, 0])
    corr = pd.read_csv(META_DIR / "rank_concordance_matrix.csv", index_col=0)
    sns.heatmap(corr, annot=True, fmt=".2f", cmap="RdYlBu_r", center=0, vmin=-1, vmax=1,
                square=True, linewidths=0.5, ax=axA, cbar_kws={"label": "Spearman rho"})
    axA.set_title("A  Cross-study signed-rank concordance", loc="left", fontweight="bold")
    plt.setp(axA.get_xticklabels(), rotation=30, ha="right", fontsize=8)
    plt.setp(axA.get_yticklabels(), fontsize=8)
    # class counts
    axB = fig.add_subplot(gs[0, 1])
    cls = pd.read_csv(META_DIR / "heterogeneity_classification.csv")
    counts = cls["class"].value_counts().reindex(
        ["I_conserved_low_het", "II_conserved_high_het", "III_context_dependent", "IV_non_reproducible"]).fillna(0)
    palette = ["#2ECC71", "#F39C12", "#3498DB", "#E74C3C"]
    axB.bar(range(len(counts)), counts.values, color=palette, ec="white")
    axB.set_xticks(range(len(counts)))
    axB.set_xticklabels(["I\nconserved\nlow-het", "II\nconserved\nhigh-het", "III\ncontext-dep", "IV\nnon-reprod"], fontsize=8)
    axB.set_ylabel("translation-factor genes")
    axB.set_title("B  Reproducibility classification", loc="left", fontweight="bold")
    for i, v in enumerate(counts.values):
        axB.text(i, v + 1, int(v), ha="center", fontsize=8)
    fig.suptitle("Figure 3. Rank concordance and gene reproducibility classes", fontsize=13, fontweight="bold")
    save(fig, "Figure3_modules")


def fig4_ml():
    lodo = pd.read_csv(ML_DIR / "leave_dataset_out_performance.csv")
    abl = pd.read_csv(ML_DIR / "ablation_results.csv")
    nc = pd.read_csv(ML_DIR / "negative_control_results.csv")
    nulld = pd.read_csv(ML_DIR / "null_distribution_random_sets.csv") if (ML_DIR / "null_distribution_random_sets.csv").exists() else None
    nsum = pd.read_csv(ML_DIR / "signature_vs_null_summary.csv") if (ML_DIR / "signature_vs_null_summary.csv").exists() else None
    perstudy = pd.read_csv(ML_DIR / "per_held_out_study_performance.csv") if (ML_DIR / "per_held_out_study_performance.csv").exists() else None
    boot = pd.read_csv(ML_DIR / "bootstrap_confidence_intervals.csv") if (ML_DIR / "bootstrap_confidence_intervals.csv").exists() else None
    ablsig = pd.read_csv(ML_DIR / "ablation_significance.csv") if (ML_DIR / "ablation_significance.csv").exists() else None

    fig = plt.figure(figsize=(15, 9.5))
    gs = fig.add_gridspec(2, 3, hspace=0.42, wspace=0.34)

    # A feature set comparison
    axA = fig.add_subplot(gs[0, 0])
    fs_order = ["all_genes", "translation_factors", "dhx_ddx_helicases", "eif_rps_rpl", "meta_signature", "random_matched"]
    best = lodo.groupby("feature_set")["balanced_accuracy"].mean().reindex(fs_order).dropna()
    bc = ["#C0392B" if f in ("translation_factors", "meta_signature") else "#4E79A7" for f in best.index]
    axA.barh([f.replace("_", " ") for f in best.index], best.values, color=bc, ec="white")
    axA.axvline(0.5, ls="--", color="#888"); axA.set_xlim(0, 1)
    axA.set_xlabel("mean LODO balanced accuracy")
    axA.set_title("A  Feature-set generalization", loc="left", fontweight="bold")
    axA.invert_yaxis()

    # B bootstrap CIs (point estimate with 95% CI)
    axB = fig.add_subplot(gs[0, 1])
    if boot is not None:
        b = boot.set_index("metric").reindex(["balanced_accuracy", "roc_auc", "pr_auc", "f1", "mcc"]).dropna()
        ys = np.arange(len(b))
        axB.errorbar(b["point_estimate"], ys,
                     xerr=[b["point_estimate"] - b["ci_lo"], b["ci_hi"] - b["point_estimate"]],
                     fmt="o", color="#2C3E50", ecolor="#E15759", elinewidth=2, capsize=4, ms=7)
        axB.set_yticks(ys); axB.set_yticklabels([m.replace("_", " ") for m in b.index], fontsize=9)
        axB.axvline(0.5, ls="--", color="#888"); axB.set_xlim(0, 1.05)
        axB.set_xlabel("score (95% bootstrap CI)")
    axB.set_title("B  Pooled LODO performance (linear SVM)", loc="left", fontweight="bold")

    # C null distribution with empirical p
    axC = fig.add_subplot(gs[0, 2])
    if nulld is not None and nsum is not None:
        vals = nulld["random_set_balanced_accuracy"].dropna()
        axC.hist(vals, bins=30, color="#BDC3C7", ec="white", alpha=0.9)
        sig = float(nsum["signature_balanced_accuracy"].iloc[0])
        axC.axvline(sig, color="#C0392B", lw=2.5, label="conserved signature")
        pct = float(nsum["percentile_vs_null"].iloc[0]); ep = float(nsum["empirical_p"].iloc[0])
        axC.text(0.04, 0.95, f"signature > {pct:.0f}% of\n1000 matched sets\nempirical p = {ep:.3g}",
                 transform=axC.transAxes, va="top", fontsize=8.5,
                 bbox=dict(boxstyle="round", fc="#FDF2E9", ec="#E67E22"))
        axC.set_xlabel("LODO balanced accuracy"); axC.set_ylabel("random matched sets")
        axC.legend(frameon=False, fontsize=8, loc="upper right")
    axC.set_title("C  Expression/variance-matched null", loc="left", fontweight="bold")

    # D ablation with significance
    axD = fig.add_subplot(gs[1, 0])
    keep = abl[abl["condition"].str.contains("top|random")].copy()
    order = ["remove_top10_signature", "remove_random10_matched", "remove_top25_signature",
             "remove_random25_matched", "remove_top50_signature", "remove_random50_matched"]
    keep = keep.set_index("condition").reindex(order).dropna().reset_index()
    colors = ["#C0392B" if "top" in c else "#BDC3C7" for c in keep["condition"]]
    axD.bar(range(len(keep)), keep["drop"], color=colors, ec="white")
    axD.set_xticks(range(len(keep)))
    axD.set_xticklabels([c.replace("remove_", "").replace("_signature", "\nsignature").replace("_matched", "\nrandom") for c in keep["condition"]], fontsize=7)
    axD.set_ylabel("drop in LODO balanced accuracy")
    if ablsig is not None:
        r25 = ablsig[ablsig["top_n_removed"] == 25]
        if not r25.empty:
            axD.text(0.5, 0.95, f"top-25 worse than {r25['signature_worse_than_pct_random'].iloc[0]:.1f}% random\n(p={r25['empirical_p'].iloc[0]:.3g})",
                     transform=axD.transAxes, va="top", ha="center", fontsize=8,
                     bbox=dict(boxstyle="round", fc="#FDEDEC", ec="#C0392B"))
    axD.set_title("D  Ablation: signature vs random matched", loc="left", fontweight="bold")

    # E negative controls
    axE = fig.add_subplot(gs[1, 1])
    nc2 = nc.copy()
    nice = {"true_labels_translation_factors": "translation factors\n(true labels)",
            "label_permutation_mean": "label permutation",
            "housekeeping_genes": "housekeeping",
            "interferon_ISG_only": "interferon (ISG)",
            "random_matched_genes_mean": "random matched",
            "study_identity_prediction_cv3_accuracy": "study-identity\n(batch probe)"}
    nc2["label"] = nc2["control"].map(nice).fillna(nc2["control"])
    nc2 = nc2.sort_values("loso_balanced_accuracy")
    bar_colors = ["#C0392B" if "true" in c else ("#7F8C8D" if "study" in c else "#4E79A7") for c in nc2["control"]]
    axE.barh(nc2["label"], nc2["loso_balanced_accuracy"], color=bar_colors, ec="white")
    axE.axvline(0.5, ls="--", color="#888"); axE.set_xlim(0, 1)
    axE.set_xlabel("LODO balanced accuracy")
    axE.set_title("E  Negative controls", loc="left", fontweight="bold")

    # F per-held-out-study performance
    axF = fig.add_subplot(gs[1, 2])
    if perstudy is not None and not perstudy.empty:
        ps = perstudy.copy()
        x = np.arange(len(ps)); w = 0.38
        axF.bar(x - w/2, ps["balanced_accuracy"], w, label="balanced acc", color="#59A14F")
        axF.bar(x + w/2, ps["roc_auc"], w, label="ROC-AUC", color="#E15759")
        axF.set_xticks(x)
        axF.set_xticklabels([f"{r.held_out_study}\n({r.virus}, n={r.n_test})" for r in ps.itertuples()], fontsize=7)
        axF.axhline(0.5, ls="--", color="#888"); axF.set_ylim(0, 1.05)
        axF.set_ylabel("score"); axF.legend(frameon=False, fontsize=8)
    axF.set_title("F  Per-held-out-study performance", loc="left", fontweight="bold")

    fig.suptitle("Figure 4. Leave-dataset-out machine learning: generalization, bootstrap CIs, matched-null significance, ablation, and controls",
                 fontsize=12.5, fontweight="bold")
    save(fig, "Figure4_ml")


def fig5_network():
    G = nx.read_graphml(NET_DIR / "poxvirus_translation_helicase_network.graphml")
    rank = pd.read_csv(SYN_DIR / "final_candidate_ranking.csv")
    fig = plt.figure(figsize=(13, 6.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.25, 1], wspace=0.22)

    axA = fig.add_subplot(gs[0, 0]); axA.set_axis_off()
    try:
        pos = nx.kamada_kawai_layout(G)
    except Exception:
        pos = nx.spring_layout(G, seed=42, k=0.9)
    node_colors = [MODULE_COLORS.get(G.nodes[n].get("module", "other"), "#95A5A6") for n in G.nodes]
    deg_c = {n: float(G.nodes[n].get("degree_centrality", 0)) for n in G.nodes}
    sizes = [260 + 2200 * deg_c[n] for n in G.nodes]
    nx.draw_networkx_edges(G, pos, alpha=0.18, ax=axA, width=0.7)
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=sizes, ax=axA, linewidths=0.6, edgecolors="white")
    # label only hub nodes (top by degree centrality) to avoid overlap
    hubs = sorted(deg_c, key=deg_c.get, reverse=True)[:12]
    nx.draw_networkx_labels(G, {n: pos[n] for n in hubs}, labels={n: n for n in hubs},
                            font_size=8, font_weight="bold", ax=axA)
    handles = [plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=c, markersize=9, label=m)
               for m, c in MODULE_COLORS.items() if m != "other"]
    axA.legend(handles=handles, frameon=False, fontsize=7.5, loc="lower left")
    axA.set_title("A  Multi-evidence translation/helicase network (hubs labeled)", loc="left", fontweight="bold")

    axB = fig.add_subplot(gs[0, 1])
    top = rank.head(15).iloc[::-1]
    bar_colors = [MODULE_COLORS.get(m, "#95A5A6") for m in top["module"]]
    axB.barh(top["gene"], top["final_evidence_score"], color=bar_colors, ec="white")
    for g in ["DHX15", "DHX29"]:
        if g in list(top["gene"]):
            i = list(top["gene"]).index(g)
            axB.get_yticklabels()[i].set_fontweight("bold")
    axB.set_xlabel("final evidence score")
    axB.set_title("B  Top multi-evidence candidates", loc="left", fontweight="bold")
    fig.suptitle("Figure 5. Network centrality and multi-evidence candidate prioritization", fontsize=13, fontweight="bold")
    save(fig, "Figure5_network")


def fig6_validation():
    rank = pd.read_csv(SYN_DIR / "final_candidate_ranking.csv")
    ext = pd.read_csv(SYN_DIR / "extended_context_effects_long.csv")
    ext["gene_symbol"] = ext["gene_symbol"].astype(str).str.upper()
    fig = plt.figure(figsize=(13, 5.4))
    gs = fig.add_gridspec(1, 3, wspace=0.34)

    # A external GSE185520 directionality for focal genes
    axA = fig.add_subplot(gs[0, 0])
    focal = ["DHX15", "DHX29", "EIF4B", "EIF3L", "DDX21", "RPL10"]
    g185 = ext[ext["context"].str.contains("GSE185520")]
    piv = g185[g185["gene_symbol"].isin(focal)].pivot_table(index="gene_symbol", columns="context", values="log2fc")
    piv = piv.reindex([f for f in focal if f in piv.index])
    if not piv.empty:
        sns.heatmap(piv, cmap="vlag", center=0, vmin=-2.5, vmax=2.5, annot=True, fmt=".1f",
                    linewidths=0.4, ax=axA, cbar_kws={"label": "log2FC"})
    axA.set_title("A  External GSE185520 validation", loc="left", fontweight="bold")
    axA.set_xlabel(""); axA.set_ylabel("")
    plt.setp(axA.get_xticklabels(), rotation=25, ha="right", fontsize=7)

    # B evidence score distribution
    axB = fig.add_subplot(gs[0, 1])
    axB.hist(rank["final_evidence_score"], bins=25, color="#4E79A7", ec="white")
    for g, c in [("DHX15", "#C0392B"), ("DHX29", "#E67E22")]:
        r = rank[rank["gene"] == g]
        if not r.empty:
            axB.axvline(r["final_evidence_score"].iloc[0], color=c, lw=2, label=g)
    axB.set_xlabel("final evidence score"); axB.set_ylabel("genes")
    axB.set_title("B  Evidence score distribution", loc="left", fontweight="bold")
    axB.legend(frameon=False, fontsize=8)

    # C conceptual model diagram
    axC = fig.add_subplot(gs[0, 2]); axC.set_axis_off()
    axC.set_xlim(0, 1); axC.set_ylim(0, 1)
    axC.set_title("C  Conceptual model", loc="left", fontweight="bold")

    def box(x, y, w, h, text, fc, fs=8):
        axC.add_patch(plt.Rectangle((x - w/2, y - h/2), w, h, fc=fc, ec="#444", lw=0.9, zorder=2))
        axC.text(x, y, text, ha="center", va="center", fontsize=fs, zorder=3)

    def arrow(x1, y1, x2, y2):
        axC.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=12,
                                      color="#555", lw=1.3, zorder=1))

    box(0.5, 0.93, 0.66, 0.1, "Poxvirus infection", "#FADBD8")
    arrow(0.5, 0.88, 0.5, 0.82)
    box(0.5, 0.77, 0.74, 0.1, "Host translational shutoff", "#FCF3CF")
    arrow(0.5, 0.72, 0.5, 0.66)
    box(0.5, 0.61, 0.86, 0.11, "eIF3 / eIF4B / 40S+60S\nribosomal proteins DOWN", "#D6EAF8")
    arrow(0.5, 0.55, 0.5, 0.49)
    box(0.5, 0.44, 0.86, 0.11, "DHX/DDX helicase\ncompensation / remodeling UP", "#FDEBD0")
    arrow(0.34, 0.385, 0.27, 0.31); arrow(0.66, 0.385, 0.73, 0.31)
    box(0.27, 0.25, 0.46, 0.12, "DHX15\nconserved axis\n(Tier 1)", "#F1948A", fs=7.5)
    box(0.73, 0.25, 0.46, 0.12, "DHX29\ncontext-dependent\nstructured-RNA (Class II)", "#AED6F1", fs=7)
    arrow(0.27, 0.19, 0.4, 0.12); arrow(0.73, 0.19, 0.6, 0.12)
    box(0.5, 0.06, 0.92, 0.09, "Experimentally testable candidates (knockdown / Ribo-seq)", "#D5F5E3", fs=7.5)

    fig.suptitle("Figure 6. External validation and integrated conceptual model", fontsize=13, fontweight="bold")
    save(fig, "Figure6_validation")


def main():
    fig1_atlas()
    fig2_meta()
    fig3_modules()
    fig4_ml()
    fig5_network()
    fig6_validation()
    print(f"Six main figures written to {FIG_DIR}")


if __name__ == "__main__":
    main()
