#!/usr/bin/env python3
"""Genome-wide, unbiased reproducibility analysis of host remodeling across poxvirus datasets.

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
META_DIR = REPO_ROOT / "results" / "meta_analysis"
OUT_DIR = REPO_ROOT / "results" / "genomewide"
FIG_DIR = OUT_DIR / "figures"

PAN_META = META_DIR / "meta_pan_poxvirus_full.csv"
VACV_META = META_DIR / "meta_vaccinia_full.csv"

# A reproducible host gene must be measured in most independent studies, agree
# in direction across them, clear genome-wide significance, and move by a
# biologically meaningful amount. These thresholds are reported, not tuned.
MIN_STUDIES = 3
FDR_CUTOFF = 0.05
EFFECT_CUTOFF = 0.5

# ---------------------------------------------------------------------------
# Curated host programs. Prefix-based families are matched by regex; pathway
# sets are explicit, literature-standard gene lists. These are transparent and
# editable; they are NOT derived from the data they are tested against.
# ---------------------------------------------------------------------------
PREFIX_PROGRAMS = {
    "RNA helicases (DEAD/DEAH-box)": r"^(DDX|DHX)\d",
    "Translation initiation (eIF)": r"^EIF",
    "Ribosomal proteins (cytoplasmic)": r"^(RPL|RPS)\d",
    "Mitochondrial ribosome": r"^(MRPL|MRPS)\d",
}

GENESET_PROGRAMS = {
    "Interferon-stimulated genes (ISG)": [
        "ISG15", "MX1", "MX2", "OAS1", "OAS2", "OAS3", "OASL", "IFIT1", "IFIT2",
        "IFIT3", "IFIT5", "IFITM1", "IFITM2", "IFITM3", "RSAD2", "IFI6", "IFI27",
        "IFI44", "IFI44L", "IFIH1", "DDX58", "STAT1", "STAT2", "IRF7", "IRF9",
        "BST2", "USP18", "HERC5", "HERC6", "XAF1", "GBP1", "GBP2", "GBP3",
        "GBP4", "GBP5", "CMPK2", "EIF2AK2", "ZBP1", "SAMD9", "SAMD9L", "MB21D1",
        "CGAS", "DDX60", "PARP9", "PARP14", "TRIM22",
    ],
    "Inflammatory / NF-kB signaling": [
        "IL6", "IL1B", "IL1A", "TNF", "CXCL8", "IL8", "CXCL1", "CXCL2", "CXCL3",
        "CXCL10", "CXCL11", "CCL2", "CCL5", "CCL20", "NFKB1", "NFKB2", "RELB",
        "REL", "NFKBIA", "NFKBIE", "NFKBIZ", "BIRC3", "TNFAIP3", "TNFAIP2",
        "TRAF1", "PTGS2", "CSF2", "CSF1", "IL23A", "TNFSF9", "TNFSF18", "ICAM1",
    ],
    "Integrated stress response / UPR": [
        "ATF3", "ATF4", "DDIT3", "CHAC1", "TRIB3", "ASNS", "PPP1R15A", "EIF2AK3",
        "EIF2AK1", "EIF2AK4", "XBP1", "HSPA5", "DDIT4", "SESN2", "STC2", "VEGFA",
        "SLC7A11", "MTHFD2", "PSAT1", "PHGDH", "CEBPB",
    ],
    "Cell cycle / mitosis": [
        "CCNB1", "CCNB2", "CCNA2", "CCNE2", "CDK1", "MKI67", "TOP2A", "BUB1",
        "BUB1B", "AURKA", "AURKB", "CDC20", "PLK1", "FOXM1", "E2F1", "CENPA",
        "CENPF", "KIF11", "TYMS", "RRM2", "PCNA",
    ],
    "Cholesterol / sterol biosynthesis": [
        "HMGCR", "HMGCS1", "LDLR", "SQLE", "DHCR7", "DHCR24", "INSIG1", "MVD",
        "FDPS", "IDI1", "MSMO1", "SC5D", "CYP51A1", "SREBF2", "FDFT1", "ACAT2",
        "LSS", "NSDHL",
    ],
    "Apoptosis / cell death": [
        "BAX", "BAK1", "BBC3", "PMAIP1", "BID", "CASP3", "CASP7", "CASP8",
        "CASP9", "TP53", "CDKN1A", "GADD45A", "FAS", "TNFRSF10A", "TNFRSF10B",
        "BCL2L11", "MCL1",
    ],
}


def load_meta(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df["gene_symbol"] = df["gene_symbol"].astype(str).str.upper()
    return df


def reproducible_set(df: pd.DataFrame) -> pd.DataFrame:
    """Genes reproducibly remodeled across independent studies (no family filter)."""
    core = df[df["k_studies"] >= MIN_STUDIES].copy()
    rep = core[
        (core["meta_FDR"] < FDR_CUTOFF)
        & (core["direction_concordance"] == 1.0)
        & (core["pooled_log2FoldChange"].abs() >= EFFECT_CUTOFF)
    ].copy()
    rep["direction"] = np.where(rep["pooled_log2FoldChange"] > 0, "up", "down")
    return rep.sort_values("pooled_log2FoldChange", key=lambda s: s.abs(), ascending=False)


def assign_programs(symbols: pd.Index) -> dict[str, set[str]]:
    """Map each curated program to the set of symbols present in `symbols`."""
    sym = pd.Index(symbols).astype(str).str.upper()
    members: dict[str, set[str]] = {}
    for name, pattern in PREFIX_PROGRAMS.items():
        rx = re.compile(pattern, re.IGNORECASE)
        members[name] = set(sym[sym.to_series().str.match(rx)])
    for name, genes in GENESET_PROGRAMS.items():
        g = {x.upper() for x in genes}
        members[name] = set(sym) & g
    return members


def hypergeom_enrichment(background: pd.Index, hits: pd.Index,
                         programs: dict[str, set[str]], label: str) -> pd.DataFrame:
    M = len(set(background))
    N = len(set(hits) & set(background))
    rows = []
    bg = set(background)
    hit = set(hits) & bg
    for name, members in programs.items():
        n = len(members & bg)
        if n == 0:
            continue
        k = len(members & hit)
        # P(X >= k) under hypergeometric null
        p = float(stats.hypergeom.sf(k - 1, M, n, N)) if k > 0 else 1.0
        expected = n * N / M if M else np.nan
        fold = (k / expected) if expected else np.nan
        rows.append({
            "program": name, "set": label, "program_size_bg": n,
            "hits_in_program": k, "expected": round(expected, 2),
            "fold_enrichment": round(fold, 2), "p_value": p,
        })
    out = pd.DataFrame(rows).sort_values("p_value")
    # BH-FDR across programs within this set
    if len(out):
        p = out["p_value"].to_numpy()
        order = np.argsort(p)
        ranked = p[order]
        m = len(ranked)
        fdr = ranked * m / (np.arange(1, m + 1))
        fdr = np.minimum.accumulate(fdr[::-1])[::-1].clip(0, 1)
        out_fdr = np.empty(m)
        out_fdr[order] = fdr
        out["FDR"] = out_fdr
    return out


def program_direction_summary(df: pd.DataFrame, programs: dict[str, set[str]],
                              min_studies: int | None = None) -> pd.DataFrame:
    """For every program, summarize reproducible direction among measured genes."""
    k_min = MIN_STUDIES if min_studies is None else min_studies
    core = df[df["k_studies"] >= k_min].copy()
    sym_to_l2fc = core.set_index("gene_symbol")["pooled_log2FoldChange"]
    sym_to_fdr = core.set_index("gene_symbol")["meta_FDR"]
    rows = []
    for name, members in programs.items():
        present = [g for g in members if g in sym_to_l2fc.index]
        if not present:
            continue
        eff = sym_to_l2fc.loc[present]
        fdr = sym_to_fdr.loc[present]
        sig = fdr < FDR_CUTOFF
        rows.append({
            "program": name,
            "genes_measured": len(present),
            "median_pooled_log2FC": round(float(eff.median()), 3),
            "mean_pooled_log2FC": round(float(eff.mean()), 3),
            "frac_up": round(float((eff > 0).mean()), 3),
            "n_sig_FDR05": int(sig.sum()),
            "n_sig_up": int(((eff > 0) & sig).sum()),
            "n_sig_down": int(((eff < 0) & sig).sum()),
        })
    return pd.DataFrame(rows).sort_values("median_pooled_log2FC")


def is_named_gene(symbols: pd.Series) -> pd.Series:
    """True for interpretable HGNC-style symbols (exclude ENSG ids, clones)."""
    s = symbols.astype(str)
    return ~s.str.match(r"^(ENSG|ENST|AC\d|AL\d|AP\d|LINC|RP\d|CTD-|CTC-|CTB-)")


def genomewide_volcano(df: pd.DataFrame, rep: pd.DataFrame, programs: dict[str, set[str]]) -> None:
    plot = df[(df["k_studies"] >= MIN_STUDIES)].dropna(
        subset=["meta_FDR", "pooled_log2FoldChange"]).copy()
    plot["mlog10"] = -np.log10(plot["meta_FDR"].clip(lower=1e-300))
    helicase = programs.get("RNA helicases (DEAD/DEAH-box)", set())

    fig, ax = plt.subplots(figsize=(9.5, 7.5))
    ns = plot[plot["meta_FDR"] >= FDR_CUTOFF]
    up = plot[(plot["meta_FDR"] < FDR_CUTOFF) & (plot["pooled_log2FoldChange"] > 0)]
    dn = plot[(plot["meta_FDR"] < FDR_CUTOFF) & (plot["pooled_log2FoldChange"] < 0)]
    ax.scatter(ns["pooled_log2FoldChange"], ns["mlog10"], s=6, c="#C7CBD1", alpha=0.35, linewidths=0)
    ax.scatter(up["pooled_log2FoldChange"], up["mlog10"], s=12, c="#C94A4A", alpha=0.55, linewidths=0, label="reproducibly up")
    ax.scatter(dn["pooled_log2FoldChange"], dn["mlog10"], s=12, c="#3E6FB6", alpha=0.55, linewidths=0, label="reproducibly down")
    hel = plot[plot["gene_symbol"].isin(helicase)]
    ax.scatter(hel["pooled_log2FoldChange"], hel["mlog10"], s=42, facecolors="none",
               edgecolors="#6B4EA0", linewidths=1.3, label="RNA helicase")
    ax.axhline(-np.log10(FDR_CUTOFF), ls="--", c="#444", lw=0.8)
    ax.axvline(0, c="#444", lw=0.8)

    # Label significant, named genes chosen by a balanced effect x significance
    # score, taken separately per direction so labels spread vertically instead
    # of piling up among the high-effect / modest-significance genes.
    named = rep[is_named_gene(rep["gene_symbol"])].merge(
        plot[["gene_symbol", "mlog10"]], on="gene_symbol", how="left")
    named = named.dropna(subset=["mlog10"])
    named["score"] = named["pooled_log2FoldChange"].abs() * named["mlog10"]
    lab = pd.concat([
        named[named["direction"] == "up"].sort_values("score", ascending=False).head(9),
        named[named["direction"] == "down"].sort_values("score", ascending=False).head(9),
    ]).drop_duplicates("gene_symbol")
    for _, r in lab.iterrows():
        ax.annotate(r["gene_symbol"], (r["pooled_log2FoldChange"], r["mlog10"]),
                    textcoords="offset points", xytext=(3, 3), fontsize=7,
                    ha="left" if r["pooled_log2FoldChange"] > 0 else "right", color="#222")
    ax.set_xlabel("pooled log2 fold change (random-effects, pan-poxvirus)")
    ax.set_ylabel("-log10 meta FDR")
    ax.set_title("Genome-wide reproducible host remodeling across independent poxvirus datasets\n"
                 "(unbiased; RNA helicases circled, not pre-selected)")
    ax.legend(frameon=False, fontsize=9, loc="upper center", ncol=2)
    fig.tight_layout()
    for suf in ("png", "pdf", "svg"):
        fig.savefig(FIG_DIR / f"genomewide_volcano.{suf}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def enrichment_barplot(enr_all: pd.DataFrame) -> None:
    d = enr_all[enr_all["set"] == "reproducible (all)"].copy()
    d = d[d["hits_in_program"] > 0].sort_values("fold_enrichment", ascending=True)
    if d.empty:
        return
    fig, ax = plt.subplots(figsize=(8.5, max(3, len(d) * 0.5)))
    colors = ["#C94A4A" if f >= 1 else "#3E6FB6" for f in d["fold_enrichment"]]
    bars = ax.barh(d["program"], d["fold_enrichment"], color=colors, alpha=0.85)
    ax.axvline(1.0, ls="--", c="#444", lw=0.8)
    for bar, (_, r) in zip(bars, d.iterrows()):
        star = "***" if r["FDR"] < 1e-3 else "**" if r["FDR"] < 1e-2 else "*" if r["FDR"] < 0.05 else ""
        ax.text(bar.get_width() + 0.05, bar.get_y() + bar.get_height() / 2,
                f"{r['hits_in_program']}/{r['program_size_bg']} {star}", va="center", fontsize=8)
    ax.set_xlabel("fold enrichment among reproducible host genes (obs/expected)")
    ax.set_title("Which host PROGRAMS are reproducibly remodeled across poxvirus datasets?\n"
                 "(* FDR<0.05, ** <0.01, *** <0.001)")
    fig.tight_layout()
    for suf in ("png", "pdf", "svg"):
        fig.savefig(FIG_DIR / f"program_enrichment.{suf}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    global MIN_STUDIES, FDR_CUTOFF, EFFECT_CUTOFF
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--min-studies", type=int, default=MIN_STUDIES)
    ap.add_argument("--fdr", type=float, default=FDR_CUTOFF)
    ap.add_argument("--effect", type=float, default=EFFECT_CUTOFF)
    args = ap.parse_args()
    MIN_STUDIES, FDR_CUTOFF, EFFECT_CUTOFF = args.min_studies, args.fdr, args.effect

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    pan = load_meta(PAN_META)
    background = pan[pan["k_studies"] >= MIN_STUDIES]["gene_symbol"]
    rep = reproducible_set(pan)
    rep.to_csv(OUT_DIR / "genomewide_reproducible_host_genes.csv", index=False)

    named = rep[is_named_gene(rep["gene_symbol"])]
    cols = ["gene_symbol", "k_studies", "pooled_log2FoldChange", "I2", "meta_FDR",
            "direction_concordance", "is_translation_factor"]
    named[named["direction"] == "up"].sort_values(
        "pooled_log2FoldChange", ascending=False).head(50)[cols].to_csv(
        OUT_DIR / "top50_reproducible_up.csv", index=False)
    named[named["direction"] == "down"].sort_values(
        "pooled_log2FoldChange").head(50)[cols].to_csv(
        OUT_DIR / "top50_reproducible_down.csv", index=False)

    programs = assign_programs(pan["gene_symbol"])

    # over-representation across programs (all / up / down)
    enr_frames = [
        hypergeom_enrichment(background, rep["gene_symbol"], programs, "reproducible (all)"),
        hypergeom_enrichment(background, rep[rep["direction"] == "up"]["gene_symbol"], programs, "reproducible up"),
        hypergeom_enrichment(background, rep[rep["direction"] == "down"]["gene_symbol"], programs, "reproducible down"),
    ]
    enr_all = pd.concat(enr_frames, ignore_index=True)
    enr_all.to_csv(OUT_DIR / "program_enrichment.csv", index=False)

    dir_summary = program_direction_summary(pan, programs)
    dir_summary.to_csv(OUT_DIR / "program_direction_summary.csv", index=False)

    # Vaccinia-infection-only robustness: the pan meta mixes true infection
    # (VacV vs mock) with Myxoma M003 *effector* contrasts, which could inflate
    # the inflammatory signal. Re-check program directions in the two Vaccinia
    # studies alone (k>=2, since only two vaccinia studies exist).
    vacv_dir = pd.DataFrame()
    if VACV_META.exists():
        vacv = load_meta(VACV_META)
        vacv_programs = assign_programs(vacv["gene_symbol"])
        vacv_dir = program_direction_summary(vacv, vacv_programs, min_studies=2)
        vacv_dir.to_csv(OUT_DIR / "program_direction_summary_vaccinia_only.csv", index=False)

    # honest contextualization: rank of headline helicases within genome-wide list
    ranked = pan[pan["k_studies"] >= MIN_STUDIES].copy()
    ranked["abs_rank_by_effect"] = ranked["pooled_log2FoldChange"].abs().rank(ascending=False).astype(int)
    ranked["rank_by_fdr"] = ranked["meta_FDR"].rank(method="min").astype(int)
    context_rows = []
    for g in ["DHX15", "DHX29", "DHX9", "DDX21", "DDX3X", "EIF4B", "EIF4E", "EIF3A"]:
        r = ranked[ranked["gene_symbol"] == g]
        if len(r):
            r = r.iloc[0]
            context_rows.append({
                "gene_symbol": g, "pooled_log2FoldChange": round(r["pooled_log2FoldChange"], 3),
                "meta_FDR": r["meta_FDR"], "abs_rank_by_effect": int(r["abs_rank_by_effect"]),
                "rank_by_fdr": int(r["rank_by_fdr"]), "of_total": int(len(ranked)),
            })
    pd.DataFrame(context_rows).to_csv(OUT_DIR / "helicase_rank_in_context.csv", index=False)

    # network-ready broad candidate list (top reproducible named genes, both directions)
    candidates = named.copy()
    candidates["evidence_rank"] = (
        candidates["pooled_log2FoldChange"].abs().rank(ascending=False)
        + (-np.log10(candidates["meta_FDR"].clip(lower=1e-300))).rank(ascending=False)
    ).rank(method="min").astype(int)
    candidates.sort_values("evidence_rank")[
        ["gene_symbol", "direction", "k_studies", "pooled_log2FoldChange",
         "I2", "meta_FDR", "is_translation_factor", "evidence_rank"]
    ].head(200).to_csv(OUT_DIR / "network_ready_candidates.csv", index=False)

    # figures
    genomewide_volcano(pan, rep, programs)
    enrichment_barplot(enr_all)

    # console report
    n_rep = len(rep)
    n_tf = int(rep["is_translation_factor"].sum())
    print("=" * 72)
    print("GENOME-WIDE REPRODUCIBILITY (pan-poxvirus, unbiased)")
    print("=" * 72)
    print(f"Background genes (>= {MIN_STUDIES} studies): {len(set(background))}")
    print(f"Reproducible host genes (FDR<{FDR_CUTOFF}, concordant, |log2FC|>={EFFECT_CUTOFF}): {n_rep}")
    print(f"  of which RNA helicases / translation factors: {n_tf} ({100*n_tf/n_rep:.1f}%)")
    print(f"  reproducibly UP: {(rep['direction']=='up').sum()} | DOWN: {(rep['direction']=='down').sum()}")
    print("\nProgram enrichment (reproducible all), top by significance:")
    show = enr_all[enr_all["set"] == "reproducible (all)"].sort_values("p_value").head(10)
    print(show[["program", "hits_in_program", "program_size_bg", "fold_enrichment", "FDR"]].to_string(index=False))
    print("\nProgram direction (median pooled log2FC; <0 = reproducibly down):")
    print(dir_summary.to_string(index=False))
    print("\nHelicase rank within the genome-wide list:")
    print(pd.DataFrame(context_rows).to_string(index=False))
    if len(vacv_dir):
        print("\nVaccinia-infection-only program direction (robustness vs effector contrasts):")
        print(vacv_dir.to_string(index=False))


if __name__ == "__main__":
    main()
