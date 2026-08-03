#!/usr/bin/env python3
"""Independent external validation of the reproducible host program using the
Matía et al. (2024, bioRxiv) Vaccinia dataset.

"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy import stats

REPO_ROOT = Path(__file__).resolve().parents[1]
EXT_DIR = REPO_ROOT / "data" / "external" / "matia2024"
META_DIR = REPO_ROOT / "results" / "meta_analysis"
GW_DIR = REPO_ROOT / "results" / "genomewide"
OUT_DIR = REPO_ROOT / "results" / "external_validation"
FIG_DIR = OUT_DIR / "figures"

SC_FILE = EXT_DIR / "matia2024_singlecell_metadata.xlsx"
PROT_FILE = EXT_DIR / "matia2024_proteomics.xlsx"

PROGRAMS = {
    "RNA helicase": (lambda g: bool(re.match(r"^(DDX|DHX)\d", g))),
    "Ribosomal protein": (lambda g: bool(re.match(r"^(RPL|RPS)\d", g))),
    "Translation initiation (eIF)": (lambda g: g.startswith("EIF")),
    "Inflammatory/NF-kB": (lambda g: g in {
        "IL6", "IL1B", "IL1A", "TNF", "CXCL8", "IL8", "CXCL1", "CXCL2", "CXCL10",
        "CCL2", "CCL5", "NFKB1", "NFKB2", "RELB", "NFKBIA", "BIRC3", "TNFAIP3",
        "PTGS2", "ICAM1"}),
    "Cell cycle/mitosis": (lambda g: g in {
        "CCNB1", "CCNB2", "CCNA2", "CDK1", "MKI67", "TOP2A", "BUB1", "AURKA",
        "AURKB", "CDC20", "PLK1", "FOXM1", "E2F1", "TYMS", "RRM2", "PCNA"}),
}


def program_of(gene: str) -> str:
    g = str(gene).upper()
    for name, fn in PROGRAMS.items():
        if fn(g):
            return name
    return "other"


# ---------------------------------------------------------------------------
# (A) Single-cell host shutoff across two cell types
# ---------------------------------------------------------------------------
def single_cell_shutoff() -> pd.DataFrame:
    rows = []
    for sheet, cell in [("metadata_HeLa", "HeLa"), ("metadata_BJ5ta", "BJ5ta")]:
        df = pd.read_excel(SC_FILE, sheet)
        state = df["infection state"].astype(str).str.lower()
        inf = df[state == "infected"]
        unin = df[state == "uninfected"]
        for metric in ["human_frac", "n_genes", "human_n_counts"]:
            if metric not in df.columns:
                continue
            a = inf[metric].dropna()
            b = unin[metric].dropna()
            if len(a) < 5 or len(b) < 5:
                continue
            u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
            rows.append({
                "cell_type": cell, "metric": metric,
                "n_infected": len(a), "n_uninfected": len(b),
                "median_infected": round(float(a.median()), 4),
                "median_uninfected": round(float(b.median()), 4),
                "log2_ratio_inf_over_uninf": round(float(np.log2((a.median() + 1e-9) / (b.median() + 1e-9))), 3),
                "mannwhitney_p": p,
            })
    out = pd.DataFrame(rows)
    out.to_csv(OUT_DIR / "matia_singlecell_shutoff_summary.csv", index=False)
    return out


def single_cell_figure() -> None:
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    for col, (sheet, cell) in enumerate([("metadata_HeLa", "HeLa"), ("metadata_BJ5ta", "BJ5ta")]):
        df = pd.read_excel(SC_FILE, sheet)
        state = df["infection state"].astype(str).str.lower()
        for row, metric, ylab in [(0, "human_frac", "host-read fraction"),
                                  (1, "n_genes", "host genes detected")]:
            ax = axes[row, col]
            data = [df.loc[state == "uninfected", metric].dropna(),
                    df.loc[state == "infected", metric].dropna()]
            bp = ax.boxplot(data, labels=["uninfected", "infected"], showfliers=False,
                            patch_artist=True, widths=0.6)
            for patch, c in zip(bp["boxes"], ["#3E6FB6", "#C94A4A"]):
                patch.set_facecolor(c); patch.set_alpha(0.7)
            ax.set_ylabel(ylab)
            if row == 0:
                ax.set_title(f"{cell} ({(state=='infected').sum()} inf / {(state=='uninfected').sum()} uninf)")
    fig.suptitle("Matía 2024: host transcriptional shutoff upon Vaccinia infection in two cell types",
                 fontsize=12)
    fig.tight_layout()
    for suf in ("png", "pdf", "svg"):
        fig.savefig(FIG_DIR / f"matia_singlecell_shutoff.{suf}", dpi=300, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# (B) Proteomics cross-modality validation
# ---------------------------------------------------------------------------
def protein_log2fc() -> pd.DataFrame:
    p = pd.read_excel(PROT_FILE, "proteomics_HeLa_host")
    vacv = ["VACV.WCE.1", "VACV.WCE.2", "VACV.WCE.3"]
    unin = ["uninf.WCE.1", "uninf.WCE.2", "uninf.WCE.3"]
    p = p.dropna(subset=vacv + unin, how="all").copy()
    # values are already log2 LFQ intensities -> log2FC is the difference of means
    p["protein_log2FC"] = p[vacv].mean(axis=1) - p[unin].mean(axis=1)
    # per-protein Welch t on the log2 intensities (need >=2 finite values per group)
    def _welch(r):
        a = r[vacv].astype(float).dropna(); b = r[unin].astype(float).dropna()
        if len(a) < 2 or len(b) < 2:
            return np.nan
        return stats.ttest_ind(a, b, equal_var=False).pvalue
    p["protein_p"] = p[vacv + unin].apply(_welch, axis=1)
    p["gene"] = p["gene"].astype(str).str.upper()
    out = p[["gene", "protein_log2FC", "protein_p"]].dropna(subset=["protein_log2FC"])
    out = out.groupby("gene", as_index=False).agg(
        protein_log2FC=("protein_log2FC", "mean"), protein_p=("protein_p", "min"))
    return out


def cross_modality(prot: pd.DataFrame) -> dict:
    # transcript reference: vaccinia-only meta (Matía is Vaccinia/HeLa -> cleanest match)
    vacv = pd.read_csv(META_DIR / "meta_vaccinia_full.csv")
    vacv["gene"] = vacv["gene_symbol"].astype(str).str.upper()
    tx = vacv[["gene", "pooled_log2FoldChange", "meta_FDR", "k_studies"]].rename(
        columns={"pooled_log2FoldChange": "transcript_log2FC"})
    merged = tx.merge(prot, on="gene", how="inner")
    merged["program"] = merged["gene"].map(program_of)
    merged.to_csv(OUT_DIR / "matia_proteomics_vs_transcriptome.csv", index=False)

    def corr(df):
        d = df.dropna(subset=["transcript_log2FC", "protein_log2FC"])
        if len(d) < 10:
            return (np.nan, np.nan, 0)
        rho, pr = stats.spearmanr(d["transcript_log2FC"], d["protein_log2FC"])
        conc = float((np.sign(d["transcript_log2FC"]) == np.sign(d["protein_log2FC"])).mean())
        return (rho, pr, conc)

    summary = {}
    rho_all, p_all, conc_all = corr(merged)
    summary["all_matched"] = {"n": int(len(merged.dropna(subset=['transcript_log2FC','protein_log2FC']))),
                              "spearman_rho": round(rho_all, 3), "spearman_p": p_all,
                              "directional_concordance": round(conc_all, 3)}
    # restrict to transcriptionally significant genes (the reproducible signal)
    sig = merged[merged["meta_FDR"] < 0.05]
    rho_s, p_s, conc_s = corr(sig)
    summary["transcript_FDR<0.05"] = {"n": int(len(sig.dropna(subset=['transcript_log2FC','protein_log2FC']))),
                                      "spearman_rho": round(rho_s, 3), "spearman_p": p_s,
                                      "directional_concordance": round(conc_s, 3)}

    # do transcript-down genes go down at protein level, and transcript-up up?
    sig_up = sig[sig["transcript_log2FC"] > 0]["protein_log2FC"].dropna()
    sig_dn = sig[sig["transcript_log2FC"] < 0]["protein_log2FC"].dropna()
    summary["protein_response_by_transcript_direction"] = {
        "transcript_up_median_protein_log2FC": round(float(sig_up.median()), 3) if len(sig_up) else None,
        "transcript_down_median_protein_log2FC": round(float(sig_dn.median()), 3) if len(sig_dn) else None,
        "mannwhitney_p": (stats.mannwhitneyu(sig_up, sig_dn).pvalue if len(sig_up) > 5 and len(sig_dn) > 5 else None),
    }

    # program-level protein response
    prog_rows = []
    for prog in [p for p in PROGRAMS] + ["other"]:
        sub = merged[merged["program"] == prog]["protein_log2FC"].dropna()
        if len(sub) >= 3:
            prog_rows.append({"program": prog, "n_proteins": len(sub),
                              "median_protein_log2FC": round(float(sub.median()), 3)})
    pd.DataFrame(prog_rows).to_csv(OUT_DIR / "matia_proteomics_program_response.csv", index=False)
    summary["program_response"] = prog_rows
    return summary, merged


def proteomics_figure(merged: pd.DataFrame) -> None:
    d = merged.dropna(subset=["transcript_log2FC", "protein_log2FC"])
    d = d[d["meta_FDR"] < 0.05]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))
    palette = {"RNA helicase": "#6B4EA0", "Ribosomal protein": "#3E6FB6",
               "Translation initiation (eIF)": "#5AA9A0", "Inflammatory/NF-kB": "#C94A4A",
               "Cell cycle/mitosis": "#E08E45", "other": "#C7CBD1"}
    for prog, sub in d.groupby("program"):
        ax1.scatter(sub["transcript_log2FC"], sub["protein_log2FC"], s=18,
                    c=palette.get(prog, "#C7CBD1"), alpha=0.6 if prog == "other" else 0.9,
                    label=prog if prog != "other" else None, linewidths=0)
    ax1.axhline(0, c="#444", lw=0.7); ax1.axvline(0, c="#444", lw=0.7)
    rho, p = stats.spearmanr(d["transcript_log2FC"], d["protein_log2FC"])
    ax1.set_xlabel("transcript pooled log2FC (Vaccinia meta-analysis)")
    ax1.set_ylabel("protein log2FC (Matía WCE, VACV/uninf)")
    ax1.set_title(f"Cross-modality concordance (transcript-significant genes)\nSpearman rho={rho:.2f}, p={p:.1e}, n={len(d)}")
    ax1.legend(frameon=False, fontsize=8, loc="upper left")

    prog = pd.read_csv(OUT_DIR / "matia_proteomics_program_response.csv")
    prog = prog[prog["program"] != "other"].sort_values("median_protein_log2FC")
    colors = ["#C94A4A" if v > 0 else "#3E6FB6" for v in prog["median_protein_log2FC"]]
    ax2.barh(prog["program"], prog["median_protein_log2FC"], color=colors, alpha=0.85)
    ax2.axvline(0, c="#444", lw=0.7)
    for y, (_, r) in enumerate(prog.iterrows()):
        ax2.text(r["median_protein_log2FC"], y, f"  n={r['n_proteins']}", va="center", fontsize=8)
    ax2.set_xlabel("median protein log2FC (VACV vs uninfected)")
    ax2.set_title("Protein-level response by host program")
    fig.tight_layout()
    for suf in ("png", "pdf", "svg"):
        fig.savefig(FIG_DIR / f"matia_cross_modality.{suf}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    argparse.ArgumentParser(description=__doc__).parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("MATÍA 2024 EXTERNAL VALIDATION (independent Vaccinia study)")
    print("=" * 70)

    sc = single_cell_shutoff()
    single_cell_figure()
    print("\n(A) Single-cell host shutoff across two cell types:")
    print(sc.to_string(index=False))

    prot = protein_log2fc()
    summary, merged = cross_modality(prot)
    proteomics_figure(merged)
    print("\n(B) Cross-modality (transcript vs protein):")
    for k, v in summary.items():
        if k != "program_response":
            print(f"  {k}: {v}")
    print("  program protein response:")
    for r in summary["program_response"]:
        print(f"    {r['program']}: median protein log2FC {r['median_protein_log2FC']} (n={r['n_proteins']})")


if __name__ == "__main__":
    main()
