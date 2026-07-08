#!/usr/bin/env python3
"""Publication-grade robustness layer for the PoxHostAtlas meta-analysis.

What this code is doing is three components on top of the random-effects meta-analysis:

1. The first thing that this is going to be doing is to do a leave-one-study-out (LOSO) meta-analysis: what this does is that is essentially
re-pool every gene after dropping each contributing study and then this will record whether the direction and the significance are
   preserved. A gene woll then be classified as being "robust" only if direction is preserved in all LOSO runs.

2. Heterogeneity classification: what I will be doing here is using Cochran's Q, I^2, tau^2, a 95% prediction
   interval, contributing-contrast count, and directional concordance, and also every
   translation/helicase gene is binned into the follwing "cans":
     Class I   - conserved, low heterogeneity
     Class II  - conserved direction, high heterogeneity
     Class III - context-dependent
     Class IV  - non-reproducible

3. Rank-based robustness: what this is doing is signed -log10(p) ranking per study, cross-study
   Spearman matrix, Kendall's W concordance, and also then a lightweight robust rank
   aggregation (RRA-style) score that goes across studies.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO_ROOT = Path(__file__).resolve().parents[1]
TABLE_DIR = REPO_ROOT / "results" / "tables"
EXPANDED_DIR = REPO_ROOT / "results" / "expanded"
META_DIR = REPO_ROOT / "results" / "meta"
OUT_DIR = REPO_ROOT / "results" / "meta_analysis"

TARGET_PATTERN = re.compile(r"^(DHX|DDX|EIF|RPS|RPL)", re.IGNORECASE)

STUDY_INPUTS = {
    "GSE278320": TABLE_DIR / "dge_results_full.csv",
    "GSE284044": EXPANDED_DIR / "GSE284044_6hpi_dge.csv",
    "GSE287860": META_DIR / "GSE287860_dge_results.csv",
    "GSE288000": META_DIR / "GSE288000_NTC_dge_results.csv",
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


def load_study(key: str) -> pd.DataFrame:
    frame = pd.read_csv(STUDY_INPUTS[key])
    frame["gene_symbol"] = frame["gene_symbol"].astype(str).str.upper()
    frame = frame[(frame["gene_symbol"].notna()) & (frame["gene_symbol"] != "NAN")]
    frame = frame[["gene_symbol", "log2FoldChange", "lfcSE", "padj", "pvalue"]].copy()
    frame = frame.dropna(subset=["log2FoldChange", "lfcSE"])
    frame = frame[frame["lfcSE"] > 0]
    frame["w"] = 1.0 / frame["lfcSE"] ** 2
    grouped = frame.groupby("gene_symbol").apply(
        lambda g: pd.Series(
            {
                "log2FoldChange": np.average(g["log2FoldChange"], weights=g["w"]),
                "lfcSE": np.sqrt(1.0 / g["w"].sum()),
                "pvalue": g["pvalue"].min(),
                "padj": g["padj"].min(),
            }
        ),
        include_groups=False,
    )
    return grouped


def dl_meta(yi: np.ndarray, vi: np.ndarray) -> dict:
    mask = np.isfinite(yi) & np.isfinite(vi) & (vi > 0)
    y, v = yi[mask], vi[mask]
    k = y.size
    if k == 0:
        return {}
    w = 1.0 / v
    sw = w.sum()
    y_fixed = float((w * y).sum() / sw)
    if k == 1:
        se = float(np.sqrt(v[0]))
        z = float(y[0] / se)
        return {
            "k": 1, "pooled": float(y[0]), "se": se, "tau2": 0.0, "Q": np.nan,
            "I2": np.nan, "z": z, "p": float(2 * stats.norm.sf(abs(z))),
            "pi_lo": np.nan, "pi_hi": np.nan, "n_up": int(y[0] > 0), "n_down": int(y[0] < 0),
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
    # 95% prediction interval (Higgins-Thompson) for k>=3
    pi_lo = pi_hi = np.nan
    if k >= 3:
        t_crit = stats.t.ppf(0.975, df)
        pi_half = t_crit * np.sqrt(tau2 + se_re ** 2)
        pi_lo, pi_hi = y_re - pi_half, y_re + pi_half
    return {
        "k": k, "pooled": y_re, "se": se_re, "tau2": tau2, "Q": Q, "I2": i2,
        "z": z, "p": float(2 * stats.norm.sf(abs(z))), "pi_lo": pi_lo, "pi_hi": pi_hi,
        "n_up": int((y > 0).sum()), "n_down": int((y < 0).sum()),
    }


def build_matrices(studies: list[str]):
    merged = None
    for key in studies:
        s = load_study(key)[["log2FoldChange", "lfcSE"]]
        s.columns = [f"y__{key}", f"se__{key}"]
        merged = s if merged is None else merged.join(s, how="outer")
    return merged


def leave_one_study_out(merged: pd.DataFrame, studies: list[str], focus_genes: list[str]) -> pd.DataFrame:
    y_cols = [f"y__{k}" for k in studies]
    se_cols = [f"se__{k}" for k in studies]
    rows = []
    sub = merged.loc[merged.index.intersection(focus_genes)]
    for gene in sub.index:
        Y = sub.loc[gene, y_cols].to_numpy(dtype="float64")
        V = sub.loc[gene, se_cols].to_numpy(dtype="float64") ** 2
        full = dl_meta(Y, V)
        if not full or full["k"] < 2:
            continue
        for i, drop in enumerate(studies):
            keep = [j for j in range(len(studies)) if j != i]
            loo = dl_meta(Y[keep], V[keep])
            if not loo:
                continue
            rows.append({
                "gene": gene,
                "full_meta_log2FC": full["pooled"],
                "full_meta_p": full["p"],
                "removed_study": drop,
                "loo_meta_log2FC": loo["pooled"],
                "loo_meta_p": loo["p"],
                "loo_k": loo["k"],
                "direction_preserved": np.sign(loo["pooled"]) == np.sign(full["pooled"]),
            })
    loo_df = pd.DataFrame(rows)
    if not loo_df.empty:
        loo_df["full_meta_FDR"] = loo_df["gene"].map(bh_adjust(loo_df.drop_duplicates("gene").set_index("gene")["full_meta_p"]))
        loo_df["significance_preserved"] = loo_df["loo_meta_p"].lt(0.05)
    loo_df.to_csv(OUT_DIR / "leave_one_study_out.csv", index=False)
    return loo_df


def classify(merged: pd.DataFrame, studies: list[str], loo_df: pd.DataFrame) -> pd.DataFrame:
    y_cols = [f"y__{k}" for k in studies]
    se_cols = [f"se__{k}" for k in studies]
    tf_index = [g for g in merged.index if TARGET_PATTERN.match(str(g))]
    robust_by_gene = {}
    if not loo_df.empty:
        robust_by_gene = loo_df.groupby("gene")["direction_preserved"].all().to_dict()
    records = []
    for gene in tf_index:
        Y = merged.loc[gene, y_cols].to_numpy(dtype="float64")
        V = merged.loc[gene, se_cols].to_numpy(dtype="float64") ** 2
        res = dl_meta(Y, V)
        if not res or res["k"] < 2:
            continue
        concordance = max(res["n_up"], res["n_down"]) / res["k"]
        i2 = res["I2"] if np.isfinite(res["I2"]) else 0.0
        loso_ok = robust_by_gene.get(gene, False)
        sig = res["p"] < 0.05
        if sig and i2 < 40 and concordance == 1.0:
            cls = "I_conserved_low_het"
        elif concordance >= 0.75 and (sig or loso_ok):
            cls = "II_conserved_high_het"
        elif concordance >= 0.5:
            cls = "III_context_dependent"
        else:
            cls = "IV_non_reproducible"
        records.append({
            "gene": gene, "k": res["k"], "pooled_log2FC": res["pooled"], "p": res["p"],
            "Q": res["Q"], "I2": i2, "tau2": res["tau2"], "pi_lo": res["pi_lo"], "pi_hi": res["pi_hi"],
            "direction_concordance": concordance, "loso_direction_robust": loso_ok,
            "class": cls,
        })
    out = pd.DataFrame(records)
    out["meta_FDR"] = bh_adjust(out.set_index("gene")["p"]).values
    out = out.sort_values(["class", "meta_FDR"])
    out.to_csv(OUT_DIR / "heterogeneity_classification.csv", index=False)
    return out


def rank_concordance(studies: list[str]) -> pd.DataFrame:
    signed = {}
    for key in studies:
        s = load_study(key)
        score = -np.log10(s["pvalue"].clip(lower=1e-300)) * np.sign(s["log2FoldChange"])
        signed[key] = score
    mat = pd.DataFrame(signed).dropna()
    corr = mat.corr(method="spearman")
    corr.to_csv(OUT_DIR / "rank_concordance_matrix.csv")
    # Kendall's W across studies on shared genes
    ranks = mat.rank(axis=0)
    n = ranks.shape[0]
    m = ranks.shape[1]
    Rj = ranks.sum(axis=1)
    S = ((Rj - Rj.mean()) ** 2).sum()
    W = 12 * S / (m ** 2 * (n ** 3 - n)) if n > 1 else np.nan
    # robust rank aggregation-style: median normalized rank of upregulated TFs
    tf = mat[mat.index.to_series().str.match(TARGET_PATTERN)]
    up_rank = tf.rank(ascending=False) / len(mat)
    rra = up_rank.median(axis=1).sort_values()
    rra.name = "median_normalized_rank"
    rra.to_csv(OUT_DIR / "robust_rank_aggregation.csv")
    pd.DataFrame([{"kendall_W": W, "n_shared_genes": n, "n_studies": m}]).to_csv(
        OUT_DIR / "rank_concordance_summary.csv", index=False
    )
    print(f"Kendall's W = {W:.3f} across {m} studies on {n} shared genes")
    return corr


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    studies = list(STUDY_INPUTS.keys())
    merged = build_matrices(studies)

    pan = pd.read_csv(OUT_DIR / "meta_pan_poxvirus_translation_factors.csv", index_col=0)
    pan.index = pan.index.astype(str).str.upper()
    focus = pan.index.tolist() + ["DHX15", "DHX29", "DDX21", "DHX9", "EIF4B", "EIF3L"]
    focus = sorted(set(focus))

    loo_df = leave_one_study_out(merged, studies, focus)
    cls = classify(merged, studies, loo_df)
    rank_concordance(studies)

    counts = cls["class"].value_counts().to_dict()
    print("Heterogeneity classes:", counts)
    for g in ["DHX15", "DHX29", "EIF4B", "DDX21"]:
        row = cls[cls["gene"] == g]
        if not row.empty:
            r = row.iloc[0]
            print(f"  {g}: {r['class']} (I2={r['I2']:.0f}%, concordance={r['direction_concordance']:.0%}, LOSO_robust={r['loso_direction_robust']})")
    print("Robust meta layer complete.")


if __name__ == "__main__":
    main()
