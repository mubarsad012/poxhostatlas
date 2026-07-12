#!/usr/bin/env python3
"""Network biology + multi-evidence candidate prioritization for PoxHostAtlas.

What this will be doing is building a host translation/helicase interaction network from STRING high-confidence
physical edges. Then what will happen is that this will compute centrality and community structure. And then it will fuse six
orthogonal evidence streams all of this into a transparent Final Evidence Score:

  meta significance + directional concordance + leave-one-study-out robustness
  + ML importance + network centrality + module membership + external validation.

The genes that are on the surface in meta-analysis AND ML AND network centrality are the
strongest and also they are the most defensible candidates.

Outputs:
  results/network/network_node_scores.csv
  results/network/poxvirus_translation_helicase_network.graphml
  results/network/community_assignments.csv
  results/synthesis/final_candidate_ranking.csv
"""

from __future__ import annotations

import re
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
META_DIR = REPO_ROOT / "results" / "meta_analysis"
MECH_DIR = REPO_ROOT / "results" / "mechanistic"
ML_DIR = REPO_ROOT / "results" / "ml"
SYN_DIR = REPO_ROOT / "results" / "synthesis"
NET_DIR = REPO_ROOT / "results" / "network"

TARGET_PATTERN = re.compile(r"^(DHX|DDX|EIF|RPS|RPL)", re.IGNORECASE)


def module_of(gene: str) -> str:
    g = gene.upper()
    if g.startswith(("DHX", "DDX")):
        return "DHX/DDX helicase"
    if g.startswith("EIF3"):
        return "eIF3 complex"
    if g.startswith("EIF4"):
        return "eIF4/initiation"
    if g.startswith("EIF"):
        return "other eIF"
    if g.startswith("RPS"):
        return "40S ribosomal"
    if g.startswith("RPL"):
        return "60S ribosomal"
    return "other"


def build_network() -> nx.Graph:
    edges = pd.read_csv(MECH_DIR / "string_network_edges.csv")
    G = nx.Graph()
    for _, r in edges.iterrows():
        a, b = str(r["preferredName_A"]).upper(), str(r["preferredName_B"]).upper()
        score = float(r["score"])
        if score >= 0.7:
            G.add_edge(a, b, weight=score)
    return G


def network_scores(G: nx.Graph) -> pd.DataFrame:
    if G.number_of_nodes() == 0:
        return pd.DataFrame()
    deg = nx.degree_centrality(G)
    btw = nx.betweenness_centrality(G, weight="weight")
    try:
        eig = nx.eigenvector_centrality(G, max_iter=1000, weight="weight")
    except nx.PowerIterationFailedConvergence:
        eig = {n: np.nan for n in G.nodes}
    communities = list(nx.algorithms.community.greedy_modularity_communities(G, weight="weight"))
    comm_map = {n: i for i, com in enumerate(communities) for n in com}
    rows = []
    for n in G.nodes:
        rows.append({
            "gene": n,
            "module": module_of(n),
            "degree_centrality": deg.get(n, 0.0),
            "betweenness_centrality": btw.get(n, 0.0),
            "eigenvector_centrality": eig.get(n, np.nan),
            "community": comm_map.get(n, -1),
            "degree": G.degree(n),
        })
    nodes = pd.DataFrame(rows).sort_values("degree_centrality", ascending=False)
    nodes.to_csv(NET_DIR / "network_node_scores.csv", index=False)
    nodes[["gene", "community", "module"]].to_csv(NET_DIR / "community_assignments.csv", index=False)
    for n in G.nodes:
        G.nodes[n]["module"] = module_of(n)
        G.nodes[n]["community"] = int(comm_map.get(n, -1))
        G.nodes[n]["degree_centrality"] = float(deg.get(n, 0.0))
    nx.write_graphml(G, NET_DIR / "poxvirus_translation_helicase_network.graphml")
    return nodes


def minmax(s: pd.Series) -> pd.Series:
    rng = s.max() - s.min()
    return (s - s.min()) / rng if rng > 0 else s * 0.0


def final_ranking(nodes: pd.DataFrame) -> pd.DataFrame:
    pan = pd.read_csv(META_DIR / "meta_pan_poxvirus_full.csv")
    pan["gene"] = pan["gene_symbol"].astype(str).str.upper()
    pan = pan[pan["gene"].str.match(TARGET_PATTERN)]
    pan = pan.set_index("gene")

    cls = pd.read_csv(META_DIR / "heterogeneity_classification.csv").set_index("gene")
    imp = pd.read_csv(ML_DIR / "feature_importance_consensus.csv").set_index("gene")
    net = nodes.set_index("gene") if not nodes.empty else pd.DataFrame()

    comp_path = SYN_DIR / "composite_evidence_score.csv"
    ext = pd.DataFrame()
    if comp_path.exists():
        ext = pd.read_csv(comp_path)
        ext["gene"] = ext["gene_symbol"].astype(str).str.upper()
        ext = ext.set_index("gene")

    genes = sorted(set(pan.index) | set(imp.index))
    rows = []
    for g in genes:
        meta = pan.loc[g] if g in pan.index else None
        meta_fdr = float(meta["meta_FDR"]) if meta is not None and pd.notna(meta["meta_FDR"]) else np.nan
        concord = float(meta["direction_concordance"]) if meta is not None else np.nan
        k = int(meta["k_studies"]) if meta is not None else 0
        loso = bool(cls.loc[g, "loso_direction_robust"]) if g in cls.index else False
        het_class = cls.loc[g, "class"] if g in cls.index else ""
        ml = float(imp.loc[g, "ml_importance_score"]) if g in imp.index else 0.0
        deg_c = float(net.loc[g, "degree_centrality"]) if (not net.empty and g in net.index) else 0.0
        ext_genes = [c for c in ext.columns if c.startswith("log2fc_GSE185520")] if not ext.empty else []
        ext_support = np.nan
        if g in getattr(ext, "index", []) and ext_genes:
            vals = ext.loc[g, ext_genes].astype(float).dropna()
            if len(vals) and meta is not None and pd.notna(meta["pooled_log2FoldChange"]):
                ext_support = float((np.sign(vals) == np.sign(meta["pooled_log2FoldChange"])).mean())
        rows.append({
            "gene": g, "module": module_of(g), "meta_FDR": meta_fdr,
            "pooled_log2FC": float(meta["pooled_log2FoldChange"]) if meta is not None else np.nan,
            "k_studies": k, "direction_concordance": concord, "het_class": het_class,
            "loso_robust": loso, "ml_importance": ml, "degree_centrality": deg_c,
            "external_support_GSE185520": ext_support,
        })
    out = pd.DataFrame(rows)

    out["s_meta"] = minmax(-np.log10(out["meta_FDR"].clip(lower=1e-300))).fillna(0)
    out["s_concord"] = out["direction_concordance"].fillna(0)
    out["s_loso"] = out["loso_robust"].astype(float)
    out["s_ml"] = minmax(out["ml_importance"]).fillna(0)
    out["s_net"] = minmax(out["degree_centrality"]).fillna(0)
    out["s_module"] = out["module"].ne("other").astype(float)
    out["s_external"] = out["external_support_GSE185520"].fillna(0)

    out["final_evidence_score"] = (
        2.5 * out["s_meta"] + 1.5 * out["s_concord"] + 1.5 * out["s_loso"]
        + 2.0 * out["s_ml"] + 1.5 * out["s_net"] + 0.5 * out["s_module"] + 1.0 * out["s_external"]
    )

    def tier(r):
        streams = sum([r["meta_FDR"] < 0.05 if pd.notna(r["meta_FDR"]) else False,
                       r["ml_importance"] > out["ml_importance"].quantile(0.75),
                       r["degree_centrality"] > 0,
                       bool(r["loso_robust"])])
        if streams >= 3:
            return "Tier 1 (triangulated)"
        if streams == 2:
            return "Tier 2 (corroborated)"
        return "Tier 3 (single-stream)"

    out["evidence_tier"] = out.apply(tier, axis=1)
    out = out.sort_values("final_evidence_score", ascending=False)
    out.to_csv(SYN_DIR / "final_candidate_ranking.csv", index=False)
    return out


def main() -> None:
    NET_DIR.mkdir(parents=True, exist_ok=True)
    SYN_DIR.mkdir(parents=True, exist_ok=True)
    G = build_network()
    print(f"Network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    nodes = network_scores(G)
    ranking = final_ranking(nodes)
    print("\nTop 12 multi-evidence candidates:")
    cols = ["gene", "module", "meta_FDR", "direction_concordance", "loso_robust",
            "ml_importance", "degree_centrality", "final_evidence_score", "evidence_tier"]
    print(ranking[cols].head(12).to_string(index=False))
    for g in ["DHX15", "DHX29"]:
        row = ranking[ranking["gene"] == g]
        if not row.empty:
            r = row.iloc[0]
            print(f"\n{g}: score={r['final_evidence_score']:.2f}, tier={r['evidence_tier']}, "
                  f"het_class={r['het_class']}, meta_FDR={r['meta_FDR']:.2g}")


if __name__ == "__main__":
    main()
