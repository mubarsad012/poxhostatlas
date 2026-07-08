#!/usr/bin/env python3
"""Genome-wide random-effects meta-analysis across independent poxvirus studies.

The cross-dataset pilot in this repository previously summarized agreement with a
simple directional vote count. This module replaces that heuristic with a formal
inverse-variance random-effects meta-analysis (DerSimonian-Laird estimator),
yielding pooled effect sizes, standard errors, z-tests, Benjamini-Hochberg FDR,
between-study heterogeneity (Cochran's Q, I^2, tau^2), and a directional
concordance audit for every gene shared across studies.

Three meta-analyses are produced:

* ``pan_poxvirus`` - four independent studies (two Vaccinia, two Myxoma-effector).
* ``vaccinia`` - the two true Vaccinia infection studies only (cleanest claim).
* ``effector`` - the two Myxoma M003 effector studies only.

Translation-factor (DHX/DDX/EIF/RPS/RPL) results are extracted and a publication
figure set is rendered (meta volcano, DHX29 forest plot across all contexts, top
conserved-factor forest panel, and the full contrast-level evidence heatmap).
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from scipy import stats

REPO_ROOT = Path(__file__).resolve().parents[1]
TABLE_DIR = REPO_ROOT / "results" / "tables"
META_DIR = REPO_ROOT / "results" / "meta"
EXPANDED_DIR = REPO_ROOT / "results" / "expanded"
OUT_DIR = REPO_ROOT / "results" / "meta_analysis"
FIG_DIR = OUT_DIR / "figures"

TARGET_PATTERN = re.compile(r"^(DHX|DDX|EIF|RPS|RPL)", re.IGNORECASE)

# Independent contrasts used as meta-analysis input units. One representative
# contrast per independent study to avoid pseudoreplication.
STUDY_INPUTS = {
    "GSE278320": {"path": TABLE_DIR / "dge_results_full.csv", "virus": "Vaccinia", "label": "GSE278320 Vaccinia (HAP1, total)"},
    "GSE284044": {"path": EXPANDED_DIR / "GSE284044_6hpi_dge.csv", "virus": "Vaccinia", "label": "GSE284044 Vaccinia (Vero, 6 hpi)"},
    "GSE287860": {"path": META_DIR / "GSE287860_dge_results.csv", "virus": "Myxoma_effector", "label": "GSE287860 Myxoma M003"},
    "GSE288000": {"path": META_DIR / "GSE288000_NTC_dge_results.csv", "virus": "Myxoma_effector", "label": "GSE288000 Myxoma M003 (NTC)"},
}

META_SETS = {
    "pan_poxvirus": ["GSE278320", "GSE284044", "GSE287860", "GSE288000"],
    "vaccinia": ["GSE278320", "GSE284044"],
    "effector": ["GSE287860", "GSE288000"],
}

# All harmonized contrasts for the full evidence heatmap / forest context.
ALL_CONTEXTS = {
    "GSE278320 VacV": TABLE_DIR / "dge_results_full.csv",
    "GSE284044 2hpi": EXPANDED_DIR / "GSE284044_2hpi_dge.csv",
    "GSE284044 6hpi": EXPANDED_DIR / "GSE284044_6hpi_dge.csv",
    "GSE284044 24hpi": EXPANDED_DIR / "GSE284044_24hpi_dge.csv",
    "GSE287860 M003": META_DIR / "GSE287860_dge_results.csv",
    "GSE288000 NTC": META_DIR / "GSE288000_NTC_dge_results.csv",
    "GSE288000 N4BP1": EXPANDED_DIR / "GSE288000_N4BP1_dge.csv",
    "GSE288000 ZC3H12A": EXPANDED_DIR / "GSE288000_ZC3H12A_dge.csv",
    "GSE288000 TRIM25": EXPANDED_DIR / "GSE288000_TRIM25_dge.csv",
    "GSE288000 dKO": EXPANDED_DIR / "GSE288000_ZC3H12A_N4BP1_dge.csv",
}


def bh_adjust(pvalues: pd.Series) -> pd.Series:
    values = pd.to_numeric(pvalues, errors="coerce")
    adjusted = pd.Series(np.nan, index=values.index, dtype="float64")
    valid = values.dropna().sort_values()
    if valid.empty:
        return adjusted
    n = len(valid)
    ranks = np.arange(1, n + 1)
    raw = valid.to_numpy() * n / ranks
    corrected = np.minimum.accumulate(raw[::-1])[::-1]
    adjusted.loc[valid.index] = np.clip(corrected, 0, 1)
    return adjusted


def load_study(symbol_key: str) -> pd.DataFrame:
    spec = STUDY_INPUTS[symbol_key]
    frame = pd.read_csv(spec["path"])
    frame["gene_symbol"] = frame["gene_symbol"].astype(str).str.upper()
    frame = frame[frame["gene_symbol"].notna() & (frame["gene_symbol"] != "NAN")]
    frame = frame[["gene_symbol", "log2FoldChange", "lfcSE", "padj", "baseMean"]].copy()
    frame = frame.dropna(subset=["log2FoldChange", "lfcSE"])
    frame = frame[frame["lfcSE"] > 0]
    # collapse duplicate symbols by inverse-variance weighting within study
    frame["w"] = 1.0 / frame["lfcSE"] ** 2
    grouped = frame.groupby("gene_symbol").apply(
        lambda g: pd.Series(
            {
                "log2FoldChange": np.average(g["log2FoldChange"], weights=g["w"]),
                "lfcSE": np.sqrt(1.0 / g["w"].sum()),
                "padj": g["padj"].min(),
                "baseMean": g["baseMean"].mean(),
            }
        ),
        include_groups=False,
    )
    grouped.columns = [f"{c}__{symbol_key}" for c in grouped.columns]
    return grouped


def random_effects_meta(yi: np.ndarray, vi: np.ndarray) -> dict:
    """DerSimonian-Laird random-effects meta-analysis for one gene.

    yi, vi are 1D arrays of effect sizes and their variances (NaN where missing).
    """
    mask = np.isfinite(yi) & np.isfinite(vi) & (vi > 0)
    y = yi[mask]
    v = vi[mask]
    k = y.size
    if k == 0:
        return {}
    w = 1.0 / v
    sw = w.sum()
    y_fixed = float((w * y).sum() / sw)
    if k == 1:
        return {
            "k_studies": 1,
            "pooled_log2FoldChange": float(y[0]),
            "pooled_SE": float(np.sqrt(v[0])),
            "tau2": 0.0,
            "Q": np.nan,
            "I2": np.nan,
            "z": float(y[0] / np.sqrt(v[0])),
            "p": float(2 * stats.norm.sf(abs(y[0] / np.sqrt(v[0])))),
            "n_up": int(y[0] > 0),
            "n_down": int(y[0] < 0),
            "direction_concordance": 1.0,
        }
    Q = float((w * (y - y_fixed) ** 2).sum())
    df = k - 1
    C = sw - (w ** 2).sum() / sw
    tau2 = max(0.0, (Q - df) / C) if C > 0 else 0.0
    w_star = 1.0 / (v + tau2)
    sw_star = w_star.sum()
    y_re = float((w_star * y).sum() / sw_star)
    se_re = float(np.sqrt(1.0 / sw_star))
    z = y_re / se_re
    i2 = max(0.0, (Q - df) / Q) * 100 if Q > 0 else 0.0
    n_up = int((y > 0).sum())
    n_down = int((y < 0).sum())
    return {
        "k_studies": k,
        "pooled_log2FoldChange": y_re,
        "pooled_SE": se_re,
        "tau2": tau2,
        "Q": Q,
        "I2": i2,
        "z": z,
        "p": float(2 * stats.norm.sf(abs(z))),
        "n_up": n_up,
        "n_down": n_down,
        "direction_concordance": max(n_up, n_down) / k,
    }


def run_meta(set_name: str, studies: list[str]) -> pd.DataFrame:
    merged = None
    for key in studies:
        s = load_study(key)
        merged = s if merged is None else merged.join(s, how="outer")
    y_cols = [f"log2FoldChange__{k}" for k in studies]
    se_cols = [f"lfcSE__{k}" for k in studies]
    Y = merged[y_cols].to_numpy(dtype="float64")
    V = merged[se_cols].to_numpy(dtype="float64") ** 2
    records = []
    for i in range(merged.shape[0]):
        res = random_effects_meta(Y[i], V[i])
        if res:
            res["gene_symbol"] = merged.index[i]
            records.append(res)
    out = pd.DataFrame.from_records(records).set_index("gene_symbol")
    # attach per-study effects for transparency
    for key in studies:
        out[f"log2fc_{key}"] = merged[f"log2FoldChange__{key}"]
        out[f"padj_{key}"] = merged[f"padj__{key}"]
    # FDR on genes observed in >=2 studies
    multi = out["k_studies"] >= 2
    out["meta_FDR"] = np.nan
    out.loc[multi, "meta_FDR"] = bh_adjust(out.loc[multi, "p"]).values
    out["is_translation_factor"] = out.index.to_series().str.match(TARGET_PATTERN)
    out = out.sort_values(["meta_FDR", "p"])
    out.to_csv(OUT_DIR / f"meta_{set_name}_full.csv")
    return out


def meta_volcano(meta: pd.DataFrame, set_name: str) -> None:
    plot = meta[meta["k_studies"] >= 3].copy() if set_name == "pan_poxvirus" else meta[meta["k_studies"] >= 2].copy()
    plot = plot.dropna(subset=["meta_FDR", "pooled_log2FoldChange"])
    plot["mlog10"] = -np.log10(plot["meta_FDR"].clip(lower=1e-300))
    plot["grp"] = "ns"
    plot.loc[(plot["meta_FDR"] < 0.05) & (plot["pooled_log2FoldChange"] > 0), "grp"] = "up"
    plot.loc[(plot["meta_FDR"] < 0.05) & (plot["pooled_log2FoldChange"] < 0), "grp"] = "down"
    plot.loc[plot["is_translation_factor"] & (plot["meta_FDR"] < 0.05), "grp"] = "TF"
    palette = {"ns": "#B6BBC2", "up": "#C94A4A", "down": "#3E6FB6", "TF": "#6B4EA0"}
    sns.set_theme(style="whitegrid", context="talk")
    fig, ax = plt.subplots(figsize=(9, 7))
    for g, fr in plot.groupby("grp"):
        ax.scatter(fr["pooled_log2FoldChange"], fr["mlog10"], s=10 if g == "ns" else 26,
                   c=palette[g], alpha=0.3 if g == "ns" else 0.85, label=g, linewidths=0)
    ax.axhline(-np.log10(0.05), ls="--", c="#333", lw=0.8)
    ax.axvline(0, c="#333", lw=0.8)
    lab = plot[plot["is_translation_factor"] & (plot["meta_FDR"] < 0.05)].sort_values("meta_FDR").head(22)
    for _, r in lab.iterrows():
        ax.text(r["pooled_log2FoldChange"], r["mlog10"], r.name, fontsize=7.5, ha="left", va="bottom")
    ax.set_xlabel("pooled log2 fold change (random effects)")
    ax.set_ylabel("-log10 meta FDR")
    ax.set_title(f"Random-effects meta-analysis volcano ({set_name.replace('_', ' ')})")
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    fig.tight_layout()
    for suf in ("png", "pdf", "svg"):
        fig.savefig(FIG_DIR / f"meta_volcano_{set_name}.{suf}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def context_effect_table() -> pd.DataFrame:
    frames = []
    for ctx, path in ALL_CONTEXTS.items():
        d = pd.read_csv(path)
        d["gene_symbol"] = d["gene_symbol"].astype(str).str.upper()
        d = d[["gene_symbol", "log2FoldChange", "lfcSE", "padj"]].copy()
        d["context"] = ctx
        frames.append(d)
    return pd.concat(frames, ignore_index=True)


def forest_plot_single(gene: str, ctx_tab: pd.DataFrame, pooled: pd.Series | None, path_stub: str, title: str) -> None:
    sub = ctx_tab[ctx_tab["gene_symbol"] == gene].copy()
    order = list(ALL_CONTEXTS.keys())
    sub["order"] = sub["context"].map({c: i for i, c in enumerate(order)})
    sub = sub.dropna(subset=["log2FoldChange", "lfcSE"]).sort_values("order", ascending=False)
    if sub.empty:
        return
    sub["lo"] = sub["log2FoldChange"] - 1.96 * sub["lfcSE"]
    sub["hi"] = sub["log2FoldChange"] + 1.96 * sub["lfcSE"]
    sns.set_theme(style="whitegrid", context="paper")
    fig, ax = plt.subplots(figsize=(7.2, max(3.2, len(sub) * 0.42 + 1)))
    ys = np.arange(len(sub))
    colors = ["#B94E48" if "Vaccinia" in c or "VacV" in c or "hpi" in c else "#4E79A7" for c in sub["context"]]
    for y, (_, r), col in zip(ys, sub.iterrows(), colors):
        ax.plot([r["lo"], r["hi"]], [y, y], color=col, lw=1.6, zorder=1)
        sig = pd.notna(r["padj"]) and r["padj"] < 0.05
        ax.scatter([r["log2FoldChange"]], [y], s=70 if sig else 45,
                   color=col, edgecolor="black" if sig else "none", linewidth=0.8, zorder=2)
    if pooled is not None and pd.notna(pooled.get("pooled_log2FoldChange", np.nan)):
        yb = -1
        lo = pooled["pooled_log2FoldChange"] - 1.96 * pooled["pooled_SE"]
        hi = pooled["pooled_log2FoldChange"] + 1.96 * pooled["pooled_SE"]
        ax.plot([lo, hi], [yb, yb], color="black", lw=2.2)
        ax.scatter([pooled["pooled_log2FoldChange"]], [yb], marker="D", s=90, color="black", zorder=3)
        ys = np.append(ys, yb)
        sub = pd.concat([sub, pd.DataFrame([{"context": "POOLED (pan-poxvirus RE)"}])], ignore_index=True)
    ax.axvline(0, color="#333", lw=0.8, ls="--")
    ax.set_yticks(sorted(ys))
    ax.set_yticklabels([c for c in sub["context"]][::-1] if pooled is None else list(sub["context"])[::-1])
    # rebuild labels robustly
    labels = list(sub["context"])
    ax.set_yticks(np.arange(-1, len(sub) - 1) if pooled is not None else np.arange(len(sub)))
    ax.set_yticklabels(labels)
    ax.set_xlabel("log2 fold change (95% CI)")
    ax.set_title(title)
    fig.tight_layout()
    for suf in ("png", "pdf", "svg"):
        fig.savefig(FIG_DIR / f"{path_stub}.{suf}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def evidence_heatmap(ctx_tab: pd.DataFrame, genes: list[str], stub: str, title: str) -> None:
    pivot = ctx_tab[ctx_tab["gene_symbol"].isin(genes)].pivot_table(
        index="gene_symbol", columns="context", values="log2FoldChange", aggfunc="mean"
    )
    pivot = pivot.reindex(columns=list(ALL_CONTEXTS.keys()))
    pivot = pivot.reindex(index=[g for g in genes if g in pivot.index])
    if pivot.empty:
        return
    sns.set_theme(style="white", context="paper")
    fig, ax = plt.subplots(figsize=(10, max(5, len(pivot) * 0.3)))
    sns.heatmap(pivot, cmap="vlag", center=0, vmin=-2.5, vmax=2.5, linewidths=0.3,
                linecolor="#E6E6E6", ax=ax, cbar_kws={"label": "log2 fold change"})
    ax.set_title(title)
    ax.set_xlabel("")
    ax.set_ylabel("")
    plt.setp(ax.get_xticklabels(), rotation=40, ha="right", fontsize=8)
    fig.tight_layout()
    for suf in ("png", "pdf", "svg"):
        fig.savefig(FIG_DIR / f"{stub}.{suf}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    metas = {name: run_meta(name, studies) for name, studies in META_SETS.items()}
    pan = metas["pan_poxvirus"]

    # translation-factor focused meta tables
    for name, meta in metas.items():
        tf = meta[meta["is_translation_factor"]].copy()
        tf.to_csv(OUT_DIR / f"meta_{name}_translation_factors.csv")

    # conserved hits summary
    summary_rows = []
    for name, meta in metas.items():
        sig = meta[(meta["meta_FDR"] < 0.05)]
        sig_tf = sig[sig["is_translation_factor"]]
        summary_rows.append({
            "meta_set": name,
            "genes_tested": int(meta["k_studies"].ge(2).sum()),
            "genes_3plus_studies": int(meta["k_studies"].ge(3).sum()),
            "sig_genes_FDR05": int(len(sig)),
            "sig_translation_factors_FDR05": int(len(sig_tf)),
            "median_I2_sig": float(sig["I2"].median()) if len(sig) else np.nan,
        })
    pd.DataFrame(summary_rows).to_csv(OUT_DIR / "meta_summary.csv", index=False)

    # figures
    for name in ("pan_poxvirus", "vaccinia"):
        meta_volcano(metas[name], name)

    ctx_tab = context_effect_table()
    ctx_tab.to_csv(OUT_DIR / "all_context_effects_long.csv", index=False)

    dhx29_pooled = pan.loc["DHX29"] if "DHX29" in pan.index else None
    forest_plot_single("DHX29", ctx_tab, dhx29_pooled, "forest_DHX29",
                       "DHX29 across poxvirus contexts (filled = FDR<0.05; red = Vaccinia)")

    # top conserved translation factors by pan meta
    top_tf = pan[pan["is_translation_factor"] & (pan["k_studies"] >= 3)].sort_values("meta_FDR").head(30)
    top_tf.to_csv(OUT_DIR / "top_conserved_translation_factors.csv")
    evidence_heatmap(ctx_tab, top_tf.index.tolist(), "evidence_heatmap_top_tf",
                     "Top conserved translation factors across all poxvirus contexts")

    # forest panels for key helicases
    for g in ["DHX29", "DHX15", "DDX21", "DHX9", "EIF4B"]:
        pooled = pan.loc[g] if g in pan.index else None
        forest_plot_single(g, ctx_tab, pooled, f"forest_{g}", f"{g}: cross-context forest plot")

    print("Meta-analysis complete.")
    print(pd.DataFrame(summary_rows).to_string(index=False))
    if dhx29_pooled is not None:
        print(f"\nDHX29 pan-poxvirus pooled log2FC={dhx29_pooled['pooled_log2FoldChange']:.3f} "
              f"FDR={dhx29_pooled['meta_FDR']:.3g} I2={dhx29_pooled['I2']:.1f}% k={int(dhx29_pooled['k_studies'])}")


if __name__ == "__main__":
    main()
