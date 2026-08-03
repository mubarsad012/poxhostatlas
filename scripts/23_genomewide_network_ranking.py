#!/usr/bin/env python3
"""Genome-wide (family-agnostic) network ranking of the most reproducibly
remodeled host genes across poxvirus datasets.

"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

REPO_ROOT = Path(__file__).resolve().parents[1]
ML_DIR = REPO_ROOT / "results" / "ml"
GW_DIR = REPO_ROOT / "results" / "genomewide"
OUT_DIR = GW_DIR / "network"
FIG_DIR = OUT_DIR / "figures"

# curated programs (mirrors scripts/22) for honest annotation of hubs
PROGRAM_PATTERNS = {
    "RNA helicase": r"^(DDX|DHX)\d",
    "Translation initiation (eIF)": r"^EIF",
    "Ribosomal protein": r"^(RPL|RPS)\d",
}
PROGRAM_SETS = {
    "Inflammatory/NF-kB": {
        "IL6", "IL1B", "IL1A", "TNF", "CXCL8", "IL8", "CXCL1", "CXCL2", "CXCL3",
        "CXCL10", "CXCL11", "CCL2", "CCL5", "CCL20", "NFKB1", "NFKB2", "RELB",
        "REL", "NFKBIA", "NFKBIE", "NFKBIZ", "BIRC3", "TNFAIP3", "TNFAIP2",
        "TRAF1", "PTGS2", "CSF2", "ICAM1", "TNFSF18", "TNFSF9",
    },
    "Interferon-stimulated": {
        "ISG15", "MX1", "MX2", "OAS1", "OAS2", "OAS3", "OASL", "IFIT1", "IFIT2",
        "IFIT3", "IFIT5", "IFITM1", "IFITM3", "RSAD2", "IFI6", "IFI27", "IFI44",
        "IFI44L", "IFIH1", "DDX58", "STAT1", "STAT2", "IRF7", "IRF9", "BST2",
        "USP18", "HERC5", "GBP1", "CMPK2", "EIF2AK2",
    },
    "Cell cycle/mitosis": {
        "CCNB1", "CCNB2", "CCNA2", "CCNE2", "CDK1", "MKI67", "TOP2A", "BUB1",
        "BUB1B", "AURKA", "AURKB", "CDC20", "PLK1", "FOXM1", "E2F1", "CENPA",
        "CENPF", "KIF11", "TYMS", "RRM2", "PCNA",
    },
}


def program_of(gene: str) -> str:
    g = gene.upper()
    for name, pat in PROGRAM_PATTERNS.items():
        if re.match(pat, g, re.IGNORECASE):
            return name
    for name, members in PROGRAM_SETS.items():
        if g in members:
            return name
    return "other"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--top-genes", type=int, default=250,
                    help="number of top reproducible host genes to include")
    ap.add_argument("--r-threshold", type=float, default=0.7,
                    help="absolute Pearson correlation to draw an edge")
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    # reproducible host genes (named, genome-wide, no family filter)
    rep = pd.read_csv(GW_DIR / "genomewide_reproducible_host_genes.csv")
    rep["gene_symbol"] = rep["gene_symbol"].astype(str).str.upper()
    named = rep[~rep["gene_symbol"].str.match(r"^(ENSG|ENST|AC\d|AL\d|AP\d|LINC|RP\d|CTD-|CTC-|CTB-)")]
    top = named.sort_values("meta_FDR").head(args.top_genes).copy()
    rep_lookup = top.set_index("gene_symbol")

    # sample-level expression (45 samples x ~21k genes, 3 studies), with study labels
    expr = pd.read_csv(ML_DIR / "sample_feature_matrix.csv", index_col=0)
    expr.columns = [c.upper() for c in expr.columns]
    meta = pd.read_csv(ML_DIR / "sample_metadata.csv")
    study_col = "Unnamed: 0" if "Unnamed: 0" in meta.columns else meta.columns[0]
    study_of = dict(zip(meta[study_col], meta["study"]))
    studies = pd.Series({s: study_of.get(s, "NA") for s in expr.index})

    genes = [g for g in top["gene_symbol"] if g in expr.columns]
    X = expr[genes].copy()

    # center within each study to remove batch baselines (keep infection covariation)
    Xc = X.copy()
    for s, idx in studies.groupby(studies).groups.items():
        rows = [i for i in idx if i in Xc.index]
        Xc.loc[rows] = Xc.loc[rows] - Xc.loc[rows].mean(axis=0)

    # correlation network across samples
    C = np.corrcoef(Xc.values.T)
    C = np.nan_to_num(C, nan=0.0)
    np.fill_diagonal(C, 0.0)

    G = nx.Graph()
    G.add_nodes_from(genes)
    n = len(genes)
    thr = args.r_threshold
    for i in range(n):
        for j in range(i + 1, n):
            r = C[i, j]
            if abs(r) >= thr:
                G.add_edge(genes[i], genes[j], weight=abs(r), sign=float(np.sign(r)))

    # centrality on the largest connected component (degree on full graph)
    deg = dict(G.degree())
    deg_c = nx.degree_centrality(G)
    if G.number_of_edges() > 0:
        comps = sorted(nx.connected_components(G), key=len, reverse=True)
        giant = G.subgraph(comps[0]).copy()
        btw = nx.betweenness_centrality(giant, weight="weight")
        try:
            eig = nx.eigenvector_centrality(giant, max_iter=2000, weight="weight")
        except nx.PowerIterationFailedConvergence:
            eig = {x: np.nan for x in giant.nodes}
        # Detect communities on POSITIVE-correlation edges only, so a module is a
        # set of genes that co-vary in the SAME direction (true co-regulation).
        # Clustering on |r| would merge strongly anti-correlated genes.
        giant_pos = nx.Graph()
        giant_pos.add_nodes_from(giant.nodes)
        for u, v, d in giant.edges(data=True):
            if d.get("sign", 1.0) > 0:
                giant_pos.add_edge(u, v, weight=d["weight"])
        if giant_pos.number_of_edges():
            communities = list(nx.algorithms.community.greedy_modularity_communities(giant_pos, weight="weight"))
            comm_map = {x: i for i, com in enumerate(communities) for x in com}
        else:
            comm_map = {}
    else:
        btw, eig, comm_map = {}, {}, {}

    rows = []
    for g in genes:
        l2fc = float(rep_lookup.loc[g, "pooled_log2FoldChange"]) if g in rep_lookup.index else np.nan
        fdr = float(rep_lookup.loc[g, "meta_FDR"]) if g in rep_lookup.index else np.nan
        rows.append({
            "gene": g,
            "program": program_of(g),
            "direction": "up" if l2fc > 0 else "down",
            "pooled_log2FC": round(l2fc, 3),
            "meta_FDR": fdr,
            "degree": deg.get(g, 0),
            "degree_centrality": round(deg_c.get(g, 0.0), 4),
            "betweenness_centrality": round(btw.get(g, 0.0), 4),
            "eigenvector_centrality": round(eig.get(g, np.nan), 4) if g in eig else np.nan,
            "community": comm_map.get(g, -1),
        })
    nodes = pd.DataFrame(rows)
    # hub score: connectivity x reproducibility (both rank-normalized)
    repro = (-np.log10(nodes["meta_FDR"].clip(lower=1e-300)))
    conn = nodes["eigenvector_centrality"].fillna(0) + nodes["degree_centrality"]
    nodes["hub_score"] = (conn.rank(pct=True) + repro.rank(pct=True)) / 2
    nodes = nodes.sort_values(["hub_score", "degree"], ascending=False)
    nodes.to_csv(OUT_DIR / "genomewide_network_hub_ranking.csv", index=False)
    nodes.head(30).to_csv(OUT_DIR / "top_host_hubs.csv", index=False)

    # community composition (which programs dominate each module)
    comm_rows = []
    for c in sorted(set(comm_map.values())):
        members = nodes[nodes["community"] == c]
        progs = members["program"].value_counts()
        comm_rows.append({
            "community": c, "n_genes": len(members),
            "dominant_program": progs.index[0] if len(progs) else "n/a",
            "frac_up": round(float((members["direction"] == "up").mean()), 2),
            "top_hubs": ", ".join(members.sort_values("degree", ascending=False)["gene"].head(6)),
        })
    pd.DataFrame(comm_rows).to_csv(OUT_DIR / "network_communities.csv", index=False)

    _draw_network(G, nodes, comm_map)
    _hub_barplot(nodes)

    n_pos = sum(1 for _, _, d in G.edges(data=True) if d.get("sign", 1.0) > 0)
    n_neg = G.number_of_edges() - n_pos
    print("=" * 70)
    print("GENOME-WIDE CO-VARIATION NETWORK (offline, family-agnostic)")
    print("=" * 70)
    print(f"nodes: {G.number_of_nodes()} | edges: {G.number_of_edges()} "
          f"({n_pos} positive / {n_neg} negative; |r|>={thr}, n=45 samples, study-centered)")
    print(f"communities (on positive co-regulation edges): {len(set(comm_map.values()))}")
    print("\nTop 20 host hubs (most central reproducibly remodeled genes, no family filter):")
    show = nodes.head(20)[["gene", "program", "direction", "degree",
                           "eigenvector_centrality", "hub_score"]]
    print(show.to_string(index=False))
    print("\nCommunity composition:")
    print(pd.DataFrame(comm_rows).to_string(index=False))


def _draw_network(G, nodes, comm_map) -> None:
    if G.number_of_edges() == 0:
        return
    comps = sorted(nx.connected_components(G), key=len, reverse=True)
    giant = G.subgraph(comps[0]).copy()
    pos = nx.spring_layout(giant, seed=42, k=0.35, weight="weight")
    ndf = nodes.set_index("gene")
    sizes = [40 + 900 * ndf.loc[x, "degree_centrality"] for x in giant.nodes]
    colors = ["#C94A4A" if ndf.loc[x, "direction"] == "up" else "#3E6FB6" for x in giant.nodes]
    fig, ax = plt.subplots(figsize=(11, 9))
    nx.draw_networkx_edges(giant, pos, alpha=0.12, width=0.6, ax=ax)
    nx.draw_networkx_nodes(giant, pos, node_size=sizes, node_color=colors, alpha=0.85,
                           linewidths=0.3, edgecolors="white", ax=ax)
    top_hubs = ndf.loc[[x for x in giant.nodes]].sort_values("degree", ascending=False).head(22).index
    nx.draw_networkx_labels(giant, pos, labels={x: x for x in top_hubs}, font_size=7.5, ax=ax)
    ax.set_title("Cross-study co-expression network of reproducibly remodeled host genes\n"
                 "(family-agnostic; red = up, blue = down; size = connectivity; top hubs labeled)")
    ax.axis("off")
    fig.tight_layout()
    for suf in ("png", "pdf", "svg"):
        fig.savefig(FIG_DIR / f"genomewide_network.{suf}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def _hub_barplot(nodes) -> None:
    d = nodes.head(20).iloc[::-1]
    fig, ax = plt.subplots(figsize=(8.5, 7))
    colors = ["#C94A4A" if x == "up" else "#3E6FB6" for x in d["direction"]]
    ax.barh(d["gene"], d["degree"], color=colors, alpha=0.85)
    for y, (_, r) in enumerate(d.iterrows()):
        ax.text(r["degree"] + 0.3, y, r["program"], va="center", fontsize=7, color="#444")
    ax.set_xlabel("network degree (co-expression connectivity)")
    ax.set_title("Top reproducibly-remodeled host hubs (red = up, blue = down)")
    fig.tight_layout()
    for suf in ("png", "pdf", "svg"):
        fig.savefig(FIG_DIR / f"top_host_hubs.{suf}", dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
