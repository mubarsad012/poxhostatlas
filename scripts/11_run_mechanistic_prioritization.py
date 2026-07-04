#!/usr/bin/env python3
"""This is for doing:; add host-curated and also mechanistic prioritization evidence for the poxvirus analysis."""

from __future__ import annotations

import argparse
import gzip
import io
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
import requests
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib import colors
from scipy.stats import hypergeom, pearsonr, spearmanr


REPO_ROOT = Path(__file__).resolve().parents[1]
PROCESSED_DIR = REPO_ROOT / "data" / "processed"
EXTERNAL_DIR = REPO_ROOT / "data" / "external"
GENCODE_DIR = EXTERNAL_DIR / "gencode"
RESULTS_DIR = REPO_ROOT / "results"
TABLE_DIR = RESULTS_DIR / "tables"
FIGURE_DIR = RESULTS_DIR / "figures"
SENSITIVITY_DIR = RESULTS_DIR / "sensitivity"
META_DIR = RESULTS_DIR / "meta"
MECH_DIR = RESULTS_DIR / "mechanistic"

GENCODE_URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz"
GENCODE_GTF = GENCODE_DIR / "gencode.v44.annotation.gtf.gz"
GENCODE_BIOTYPES = PROCESSED_DIR / "gencode_v44_gene_biotypes.csv"
STRING_API = "https://string-db.org/api/tsv/network"
TARGET_PATTERN = re.compile(r"^(DHX|DDX|EIF|RPS|RPL)", re.IGNORECASE)
PSEUDOGENE_SYMBOL_PATTERN = re.compile(r"^(RPL|RPS|EIF|DHX|DDX).*(?:P\d+|AP\d+)$", re.IGNORECASE)
PSEUDOGENE_SYMBOL_EXCEPTIONS = {
    "EIF4EBP1",
    "EIF4EBP2",
    "EIF4EBP3",
    "RPLP0",
    "RPLP1",
    "RPLP2",
}


CURATED_MODULES: dict[str, set[str]] = {
    "DHX_DDX_RNA_helicases": {
        "DHX9",
        "DHX15",
        "DHX29",
        "DHX30",
        "DHX33",
        "DHX36",
        "DHX37",
        "DHX57",
        "DDX3X",
        "DDX5",
        "DDX17",
        "DDX21",
        "DDX27",
        "DDX39A",
        "DDX39B",
        "DDX41",
        "DDX46",
        "DDX54",
        "DDX56",
    },
    "EIF3_complex": {f"EIF3{letter}" for letter in "ABCDEFGHIJKLM"},
    "EIF2_initiator_axis": {
        "EIF1",
        "EIF1AX",
        "EIF2A",
        "EIF2S1",
        "EIF2S2",
        "EIF2S3",
        "EIF2B1",
        "EIF2B2",
        "EIF2B3",
        "EIF2B4",
        "EIF2B5",
        "EIF5",
        "EIF5B",
    },
    "EIF4_cap_scanning_axis": {
        "EIF4A1",
        "EIF4A2",
        "EIF4B",
        "EIF4E",
        "EIF4E2",
        "EIF4EBP1",
        "EIF4EBP2",
        "EIF4EBP3",
        "EIF4ENIF1",
        "EIF4G1",
        "EIF4G2",
        "EIF4H",
    },
    "40S_ribosomal_subunits": {
        "RPS2",
        "RPS3",
        "RPS3A",
        "RPS4X",
        "RPS5",
        "RPS6",
        "RPS7",
        "RPS8",
        "RPS9",
        "RPS10",
        "RPS11",
        "RPS12",
        "RPS13",
        "RPS14",
        "RPS15",
        "RPS15A",
        "RPS16",
        "RPS17",
        "RPS18",
        "RPS19",
        "RPS20",
        "RPS21",
        "RPS23",
        "RPS24",
        "RPS25",
        "RPS26",
        "RPS27",
        "RPS27A",
        "RPS28",
        "RPS29",
        "RPSA",
    },
    "60S_ribosomal_subunits": {
        "RPL3",
        "RPL4",
        "RPL5",
        "RPL6",
        "RPL7",
        "RPL7A",
        "RPL8",
        "RPL9",
        "RPL10",
        "RPL10A",
        "RPL11",
        "RPL12",
        "RPL13",
        "RPL13A",
        "RPL14",
        "RPL15",
        "RPL17",
        "RPL18",
        "RPL18A",
        "RPL19",
        "RPL21",
        "RPL22",
        "RPL23",
        "RPL23A",
        "RPL24",
        "RPL26",
        "RPL27",
        "RPL27A",
        "RPL28",
        "RPL29",
        "RPL30",
        "RPL31",
        "RPL32",
        "RPL34",
        "RPL35",
        "RPL35A",
        "RPL36",
        "RPL36A",
        "RPL37",
        "RPL37A",
        "RPL38",
        "RPL39",
        "RPL40",
        "RPLP0",
        "RPLP1",
        "RPLP2",
    },
}
CURATED_MODULES["Translation_initiation_broad"] = (
    CURATED_MODULES["DHX_DDX_RNA_helicases"]
    | CURATED_MODULES["EIF3_complex"]
    | CURATED_MODULES["EIF2_initiator_axis"]
    | CURATED_MODULES["EIF4_cap_scanning_axis"]
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--padj", type=float, default=0.05)
    parser.add_argument("--effect-size", type=float, default=0.5)
    parser.add_argument("--top-n", type=int, default=30)
    parser.add_argument("--force-gencode", action="store_true")
    parser.add_argument("--force-string", action="store_true")
    parser.add_argument("--skip-string", action="store_true")
    return parser.parse_args()


def download_file(url: str, destination: Path, force: bool = False) -> None:
    if destination.exists() and destination.stat().st_size > 0 and not force:
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    tmp = destination.with_suffix(destination.suffix + ".tmp")
    with requests.get(url, stream=True, timeout=180) as response:
        response.raise_for_status()
        with tmp.open("wb") as handle:
            for chunk in response.iter_content(1024 * 1024):
                if chunk:
                    handle.write(chunk)
    tmp.replace(destination)


def parse_gtf_attributes(attribute_text: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for item in attribute_text.strip().split(";"):
        item = item.strip()
        if not item:
            continue
        if " " not in item:
            continue
        key, value = item.split(" ", 1)
        attributes[key] = value.strip().strip('"')
    return attributes


def build_gencode_biotypes(force: bool = False) -> pd.DataFrame:
    if GENCODE_BIOTYPES.exists() and GENCODE_BIOTYPES.stat().st_size > 0 and not force:
        return pd.read_csv(GENCODE_BIOTYPES, dtype="string")

    download_file(GENCODE_URL, GENCODE_GTF, force=force)
    rows = []
    with gzip.open(GENCODE_GTF, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            attrs = parse_gtf_attributes(fields[8])
            gene_id = attrs.get("gene_id")
            if not gene_id:
                continue
            rows.append(
                {
                    "gene_id": gene_id,
                    "gene_id_base": gene_id.split(".", 1)[0],
                    "gencode_symbol": attrs.get("gene_name", ""),
                    "gene_type": attrs.get("gene_type", attrs.get("gene_biotype", "")),
                    "gencode_source": fields[1],
                    "gencode_release": "v44",
                    "gencode_url": GENCODE_URL,
                }
            )
    biotypes = pd.DataFrame(rows).drop_duplicates("gene_id_base")
    GENCODE_BIOTYPES.parent.mkdir(parents=True, exist_ok=True)
    biotypes.to_csv(GENCODE_BIOTYPES, index=False)
    return biotypes


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


def family_from_symbol(symbol: str) -> str:
    for prefix in ("DHX", "DDX", "EIF", "RPS", "RPL"):
        if symbol.upper().startswith(prefix):
            return prefix
    return "other"


def pseudogene_symbol_flag(symbol: str) -> bool:
    symbol = str(symbol)
    if symbol in PSEUDOGENE_SYMBOL_EXCEPTIONS:
        return False
    return bool(PSEUDOGENE_SYMBOL_PATTERN.match(symbol))


def load_annotated_results(padj: float, effect_size: float, force_gencode: bool) -> pd.DataFrame:
    results = pd.read_csv(TABLE_DIR / "dge_results_full.csv")
    biotypes = build_gencode_biotypes(force=force_gencode)
    results["gene_id_base"] = results["gene_id"].astype(str).str.split(".", n=1).str[0]
    annotated = results.merge(
        biotypes[["gene_id_base", "gencode_symbol", "gene_type", "gencode_release"]],
        on="gene_id_base",
        how="left",
    )
    annotated["host_gene"] = ~annotated["gene_id"].astype(str).str.startswith("gene-")
    annotated["target_family"] = annotated["gene_symbol"].fillna("").map(family_from_symbol)
    annotated["translation_target"] = annotated["gene_symbol"].fillna("").str.match(TARGET_PATTERN)
    annotated["protein_coding"] = annotated["gene_type"].fillna("").eq("protein_coding")
    annotated["gencode_pseudogene"] = annotated["gene_type"].fillna("").str.contains("pseudogene", case=False)
    annotated["symbol_pseudogene_like"] = annotated["gene_symbol"].fillna("").map(pseudogene_symbol_flag)
    annotated["curated_host_protein_coding_target"] = (
        annotated["host_gene"]
        & annotated["translation_target"]
        & annotated["protein_coding"]
        & ~annotated["gencode_pseudogene"]
        & ~annotated["symbol_pseudogene_like"]
    )
    annotated["significant"] = annotated["padj"].notna() & (annotated["padj"] < padj)
    annotated["passes_effect_size"] = annotated["log2FoldChange"].abs() >= effect_size
    annotated["significant_with_effect"] = annotated["significant"] & annotated["passes_effect_size"]
    return annotated


def viral_load_table() -> pd.DataFrame:
    normalized = pd.read_csv(PROCESSED_DIR / "normalized_counts.csv").set_index("sample_id")
    metadata = pd.read_csv(PROCESSED_DIR / "model_ready_metadata.csv")
    viral_cols = [col for col in normalized.columns if col.startswith("gene-VAC")]
    if not viral_cols:
        raise SystemExit("No Vaccinia viral genes found in normalized count matrix.")
    load = pd.DataFrame(
        {
            "sample_id": normalized.index,
            "viral_gene_count": len(viral_cols),
            "viral_normalized_sum": normalized[viral_cols].sum(axis=1).to_numpy(),
            "viral_normalized_mean": normalized[viral_cols].mean(axis=1).to_numpy(),
        }
    )
    load["log10_viral_normalized_sum_plus1"] = np.log10(load["viral_normalized_sum"] + 1)
    return load.merge(metadata, on="sample_id", how="left")


def correlate_host_targets_with_viral_load(annotated: pd.DataFrame) -> pd.DataFrame:
    normalized = pd.read_csv(PROCESSED_DIR / "normalized_counts.csv").set_index("sample_id")
    load = viral_load_table().set_index("sample_id")
    infected_samples = load[load["infection"].eq("VacV")].index.intersection(normalized.index)
    x = load.loc[infected_samples, "log10_viral_normalized_sum_plus1"].astype(float)
    rows = []
    target_rows = annotated[annotated["curated_host_protein_coding_target"]].copy()
    for _, row in target_rows.iterrows():
        gene_id = row["gene_id"]
        if gene_id not in normalized.columns:
            continue
        y = np.log2(normalized.loc[infected_samples, gene_id].astype(float) + 1)
        if x.nunique() < 2 or y.nunique() < 2 or len(y.dropna()) < 3:
            pearson_r = pearson_p = spearman_r = spearman_p = np.nan
        else:
            pearson_r, pearson_p = pearsonr(x, y)
            spearman_r, spearman_p = spearmanr(x, y)
        rows.append(
            {
                "gene_id": gene_id,
                "gene_symbol": row["gene_symbol"],
                "target_family": row["target_family"],
                "n_infected_samples": len(infected_samples),
                "pearson_r_viral_load": pearson_r,
                "pearson_p_viral_load": pearson_p,
                "spearman_r_viral_load": spearman_r,
                "spearman_p_viral_load": spearman_p,
            }
        )
    frame = pd.DataFrame(rows)
    if not frame.empty:
        frame["viral_load_correlation_q"] = bh_adjust(frame["pearson_p_viral_load"])
        frame = frame.sort_values(["viral_load_correlation_q", "pearson_p_viral_load", "gene_symbol"])
    return frame


def module_enrichment(annotated: pd.DataFrame, padj: float, effect_size: float) -> pd.DataFrame:
    universe = set(
        annotated.loc[annotated["host_gene"] & annotated["protein_coding"], "gene_symbol"].dropna().astype(str)
    )
    queries = {
        "all_significant_host_protein_coding": set(
            annotated.loc[annotated["host_gene"] & annotated["protein_coding"] & annotated["significant"], "gene_symbol"]
        ),
        "upregulated_significant_host_protein_coding": set(
            annotated.loc[
                annotated["host_gene"]
                & annotated["protein_coding"]
                & annotated["significant"]
                & (annotated["log2FoldChange"] > 0),
                "gene_symbol",
            ]
        ),
        "downregulated_significant_host_protein_coding": set(
            annotated.loc[
                annotated["host_gene"]
                & annotated["protein_coding"]
                & annotated["significant"]
                & (annotated["log2FoldChange"] < 0),
                "gene_symbol",
            ]
        ),
        "translation_targets_significant_with_effect": set(
            annotated.loc[
                annotated["curated_host_protein_coding_target"]
                & annotated["significant"]
                & (annotated["log2FoldChange"].abs() >= effect_size),
                "gene_symbol",
            ]
        ),
    }
    rows = []
    n_universe = len(universe)
    for query_name, query in queries.items():
        query = query & universe
        for module_name, module_genes in CURATED_MODULES.items():
            module = module_genes & universe
            overlap = sorted(query & module)
            pvalue = hypergeom.sf(len(overlap) - 1, n_universe, len(module), len(query)) if query and module else np.nan
            rows.append(
                {
                    "query_set": query_name,
                    "module": module_name,
                    "universe_size": n_universe,
                    "query_size": len(query),
                    "module_size": len(module),
                    "overlap_size": len(overlap),
                    "pvalue": pvalue,
                    "overlap_genes": ";".join(overlap),
                }
            )
    frame = pd.DataFrame(rows)
    frame["padj"] = bh_adjust(frame["pvalue"])
    return frame.sort_values(["padj", "pvalue", "query_set", "module"])


def build_priority_table(
    annotated: pd.DataFrame,
    correlations: pd.DataFrame,
    padj: float,
    effect_size: float,
) -> pd.DataFrame:
    curated = annotated[annotated["curated_host_protein_coding_target"]].copy()
    comparison_path = SENSITIVITY_DIR / "sensitivity_primary_vs_strict_translation.csv"
    if comparison_path.exists():
        comparison = pd.read_csv(comparison_path)
        curated = curated.merge(
            comparison[
                [
                    "gene_symbol",
                    "log2FoldChange_strict",
                    "padj_strict",
                    "significant_strict",
                    "direction_agreement",
                ]
            ].drop_duplicates("gene_symbol"),
            on="gene_symbol",
            how="left",
        )
    meta_path = META_DIR / "poxvirus_translation_factor_meta_signature.csv"
    if meta_path.exists():
        meta = pd.read_csv(meta_path)
        curated = curated.merge(meta, on="gene_symbol", how="left", suffixes=("", "_meta"))
    if not correlations.empty:
        curated = curated.merge(correlations, on=["gene_id", "gene_symbol", "target_family"], how="left")

    for column, default in {
        "significant_dataset_count": 0,
        "dataset_count": 1,
        "direction_consistency": 0,
        "significant_strict": False,
        "direction_agreement": False,
        "pearson_r_viral_load": np.nan,
    }.items():
        if column not in curated.columns:
            curated[column] = default
        if isinstance(default, bool):
            curated[column] = curated[column].where(curated[column].notna(), default).astype(bool)
        else:
            curated[column] = curated[column].fillna(default)

    effect_component = np.minimum(curated["log2FoldChange"].abs(), 4) / 4
    viral_component = curated["pearson_r_viral_load"].abs().fillna(0)
    base_mean_component = np.log10(curated["baseMean"].fillna(0) + 1) / np.log10(curated["baseMean"].fillna(0).max() + 1)
    curated["candidate_priority_score"] = (
        curated["significant"].astype(float) * 2.0
        + curated["passes_effect_size"].astype(float) * 0.75
        + curated["significant_strict"].astype(float) * 1.0
        + curated["direction_agreement"].astype(float) * 0.5
        + curated["significant_dataset_count"].astype(float) * 1.0
        + curated["direction_consistency"].astype(float) * 0.75
        + effect_component
        + viral_component * 0.75
        + base_mean_component * 0.25
    )
    curated["priority_interpretation"] = np.select(
        [
            curated["gene_symbol"].eq("DHX29"),
            curated["significant_dataset_count"].ge(3),
            curated["significant"] & curated["significant_strict"] & curated["passes_effect_size"],
            curated["significant"],
        ],
        [
            "mentor-focal DHX29 candidate",
            "cross-dataset conserved translation factor",
            "strict-sensitivity-supported primary target",
            "primary significant translation target",
        ],
        default="audited host protein-coding translation target",
    )
    return curated.sort_values(
        ["candidate_priority_score", "significant_dataset_count", "padj"],
        ascending=[False, False, True],
    )


def fetch_string_network(symbols: list[str], force: bool = False, skip: bool = False) -> pd.DataFrame:
    edge_path = MECH_DIR / "string_network_edges.csv"
    if skip:
        return pd.DataFrame()
    if edge_path.exists() and edge_path.stat().st_size > 0 and not force:
        return pd.read_csv(edge_path)
    if not symbols:
        return pd.DataFrame()
    payload = {
        "identifiers": "\r".join(symbols),
        "species": 9606,
        "required_score": 700,
        "network_type": "physical",
        "caller_identity": "openbio_sri_poxvirus_transcriptomics",
    }
    try:
        response = requests.post(STRING_API, data=payload, timeout=120)
        response.raise_for_status()
        frame = pd.read_csv(io.StringIO(response.text), sep="\t")
    except Exception as exc:
        frame = pd.DataFrame(
            [
                {
                    "string_api_status": "failed",
                    "string_api_error": str(exc),
                    "input_symbols": ";".join(symbols),
                    "string_api_url": STRING_API,
                }
            ]
        )
    frame.to_csv(edge_path, index=False)
    return frame


def build_string_nodes(edges: pd.DataFrame, priority: pd.DataFrame) -> pd.DataFrame:
    if edges.empty or "preferredName_A" not in edges.columns:
        return pd.DataFrame()
    names = sorted(set(edges["preferredName_A"].astype(str)) | set(edges["preferredName_B"].astype(str)))
    lookup = priority.drop_duplicates("gene_symbol").set_index("gene_symbol")
    rows = []
    for name in names:
        row = lookup.loc[name] if name in lookup.index else pd.Series(dtype="object")
        degree = int((edges["preferredName_A"].eq(name) | edges["preferredName_B"].eq(name)).sum())
        rows.append(
            {
                "gene_symbol": name,
                "degree_in_string_subnetwork": degree,
                "log2FoldChange": row.get("log2FoldChange", np.nan),
                "padj": row.get("padj", np.nan),
                "candidate_priority_score": row.get("candidate_priority_score", np.nan),
                "target_family": row.get("target_family", "external_or_unranked"),
            }
        )
    nodes = pd.DataFrame(rows).sort_values(["degree_in_string_subnetwork", "candidate_priority_score"], ascending=False)
    nodes.to_csv(MECH_DIR / "string_network_nodes.csv", index=False)
    return nodes


def write_summary_metrics(annotated: pd.DataFrame, priority: pd.DataFrame, enrichment: pd.DataFrame) -> None:
    dhx29 = priority[priority["gene_symbol"].eq("DHX29")]
    helicases = priority[priority["target_family"].isin(["DHX", "DDX"])].reset_index(drop=True)
    dhx29_helicase_rank = "absent"
    if not dhx29.empty and not helicases.empty:
        helicase_match = helicases[helicases["gene_symbol"].eq("DHX29")]
        if not helicase_match.empty:
            dhx29_helicase_rank = int(helicase_match.index[0] + 1)
    rows = [
        {"metric": "host_model_genes", "value": int(annotated["host_gene"].sum())},
        {"metric": "viral_model_genes", "value": int((~annotated["host_gene"]).sum())},
        {"metric": "host_protein_coding_model_genes", "value": int((annotated["host_gene"] & annotated["protein_coding"]).sum())},
        {"metric": "curated_host_protein_coding_translation_targets", "value": int(len(priority))},
        {
            "metric": "curated_significant_with_abs_log2fc_ge_threshold",
            "value": int((priority["significant"] & priority["passes_effect_size"]).sum()),
        },
        {"metric": "dhx29_priority_rank", "value": int(dhx29.index[0] + 1) if not dhx29.empty else "absent"},
        {"metric": "dhx29_helicase_priority_rank", "value": dhx29_helicase_rank},
        {
            "metric": "module_enrichment_min_padj",
            "value": float(enrichment["padj"].min()) if not enrichment.empty else np.nan,
        },
    ]
    pd.DataFrame(rows).to_csv(MECH_DIR / "mechanistic_summary_metrics.csv", index=False)


def save_host_curated_volcano(annotated: pd.DataFrame, priority: pd.DataFrame, padj: float, top_n: int) -> None:
    host = annotated[annotated["host_gene"]].copy()
    host["minus_log10_padj"] = -np.log10(host["padj"].clip(lower=1e-300))
    host["plot_group"] = "not significant"
    host.loc[host["significant"] & (host["log2FoldChange"] > 0), "plot_group"] = "upregulated"
    host.loc[host["significant"] & (host["log2FoldChange"] < 0), "plot_group"] = "downregulated"
    curated_symbols = set(priority.loc[priority["significant"] & priority["passes_effect_size"], "gene_symbol"])
    host.loc[host["gene_symbol"].isin(curated_symbols), "plot_group"] = "curated translation target"
    palette = {
        "not significant": "#A7ADB4",
        "upregulated": "#C94A4A",
        "downregulated": "#3E6FB6",
        "curated translation target": "#6B4EA0",
    }
    sns.set_theme(style="whitegrid", context="talk")
    fig, ax = plt.subplots(figsize=(10, 7))
    for group, frame in host.groupby("plot_group"):
        ax.scatter(
            frame["log2FoldChange"],
            frame["minus_log10_padj"],
            s=14 if group == "not significant" else 22,
            alpha=0.35 if group == "not significant" else 0.85,
            c=palette[group],
            label=group,
            linewidths=0,
        )
    ax.axvline(0, color="#333333", lw=0.8)
    ax.axhline(-math.log10(padj), color="#333333", lw=0.8, ls="--")
    labels = ["DHX29", "DHX15", "DHX57", "EIF3F", "EIF3L", "EIF4B", "EIF5", "EIF1AD"]
    labels.extend(priority.head(top_n)["gene_symbol"].tolist())
    for _, row in host[host["gene_symbol"].isin(dict.fromkeys(labels))].iterrows():
        ax.text(row["log2FoldChange"], row["minus_log10_padj"], row["gene_symbol"], fontsize=8, ha="left", va="bottom")
    ax.set_title("Host-only curated translation-factor volcano")
    ax.set_xlabel("log2 fold change (VacV vs mock)")
    ax.set_ylabel("-log10 adjusted p-value")
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    fig.tight_layout()
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(MECH_DIR / f"host_only_curated_volcano.{suffix}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_dhx29_viral_load_plot(priority: pd.DataFrame) -> None:
    normalized = pd.read_csv(PROCESSED_DIR / "normalized_counts.csv").set_index("sample_id")
    load = viral_load_table().set_index("sample_id")
    dhx29 = priority[priority["gene_symbol"].eq("DHX29")]
    if dhx29.empty:
        return
    gene_id = dhx29.iloc[0]["gene_id"]
    if gene_id not in normalized.columns:
        return
    plot = load.copy()
    plot["DHX29_log2_normalized_count_plus1"] = np.log2(normalized[gene_id].astype(float) + 1)
    infected = plot[plot["infection"].eq("VacV")]
    r, p = pearsonr(
        infected["log10_viral_normalized_sum_plus1"],
        infected["DHX29_log2_normalized_count_plus1"],
    )
    sns.set_theme(style="whitegrid", context="talk")
    fig, ax = plt.subplots(figsize=(7, 5))
    sns.scatterplot(
        data=plot.reset_index(),
        x="log10_viral_normalized_sum_plus1",
        y="DHX29_log2_normalized_count_plus1",
        hue="infection",
        style="infection",
        s=90,
        ax=ax,
        palette={"mock": "#7A7A7A", "VacV": "#B94E48"},
    )
    if len(infected) >= 3:
        sns.regplot(
            data=infected,
            x="log10_viral_normalized_sum_plus1",
            y="DHX29_log2_normalized_count_plus1",
            scatter=False,
            ax=ax,
            color="#B94E48",
        )
    ax.set_title("DHX29 transcript abundance versus viral RNA burden")
    ax.set_xlabel("log10 summed normalized Vaccinia RNA + 1")
    ax.set_ylabel("log2 normalized DHX29 count + 1")
    ax.text(0.04, 0.96, f"VacV-only Pearson r={r:.2f}, p={p:.3g}", transform=ax.transAxes, va="top", fontsize=10)
    fig.tight_layout()
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(MECH_DIR / f"dhx29_viral_load_correlation.{suffix}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_module_enrichment_plot(enrichment: pd.DataFrame) -> None:
    plot = enrichment[enrichment["overlap_size"] > 0].copy().sort_values("padj").head(12)
    if plot.empty:
        return
    plot["minus_log10_padj"] = -np.log10(plot["padj"].clip(lower=1e-300))
    plot["label"] = plot["query_set"].str.replace("_", " ") + " | " + plot["module"].str.replace("_", " ")
    sns.set_theme(style="whitegrid", context="paper")
    fig, ax = plt.subplots(figsize=(9, max(4, len(plot) * 0.35)))
    sns.barplot(data=plot, x="minus_log10_padj", y="label", hue="overlap_size", dodge=False, palette="viridis", ax=ax)
    ax.set_xlabel("-log10 adjusted enrichment p-value")
    ax.set_ylabel("")
    ax.set_title("Curated translation-module enrichment")
    ax.legend(title="Overlap genes", frameon=False)
    fig.tight_layout()
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(MECH_DIR / f"module_enrichment_barplot.{suffix}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_candidate_evidence_matrix(priority: pd.DataFrame, top_n: int) -> None:
    columns = [
        "log2FoldChange",
        "log2FoldChange_strict",
        "log2fc_GSE287860",
        "log2fc_GSE288000_NTC",
        "pearson_r_viral_load",
    ]
    available = [col for col in columns if col in priority.columns]
    plot = priority.head(top_n).set_index("gene_symbol")[available].apply(pd.to_numeric, errors="coerce")
    plot = plot.rename(
        columns={
            "log2FoldChange": "GSE278320 primary log2FC",
            "log2FoldChange_strict": "strict-title log2FC",
            "log2fc_GSE287860": "GSE287860 log2FC",
            "log2fc_GSE288000_NTC": "GSE288000_NTC log2FC",
            "pearson_r_viral_load": "viral-load r",
        }
    )
    sns.set_theme(style="white", context="paper")
    fig, ax = plt.subplots(figsize=(9, max(5, len(plot) * 0.28)))
    sns.heatmap(plot, cmap="vlag", center=0, linewidths=0.2, linecolor="#DDDDDD", ax=ax, cbar_kws={"label": "standardized evidence value"})
    ax.set_title("Top mechanistic candidates across evidence layers")
    ax.set_xlabel("")
    ax.set_ylabel("")
    fig.tight_layout()
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(MECH_DIR / f"candidate_evidence_matrix.{suffix}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_helicase_priority_plot(priority: pd.DataFrame) -> None:
    helicases = priority[priority["target_family"].isin(["DHX", "DDX"])].head(25).copy()
    if helicases.empty:
        return
    helicases["is_dhx29"] = helicases["gene_symbol"].eq("DHX29")
    sns.set_theme(style="whitegrid", context="paper")
    fig, ax = plt.subplots(figsize=(8, max(4, len(helicases) * 0.28)))
    bar_colors = np.where(helicases["is_dhx29"], "#B94E48", "#4E79A7")
    ax.barh(helicases["gene_symbol"], helicases["candidate_priority_score"], color=bar_colors)
    ax.invert_yaxis()
    ax.set_xlabel("mechanistic priority score")
    ax.set_ylabel("")
    ax.set_title("DHX/DDX helicase prioritization")
    for _, row in helicases.iterrows():
        ax.text(
            row["candidate_priority_score"] + 0.05,
            row["gene_symbol"],
            f"log2FC {row['log2FoldChange']:.2f}; meta sig {int(row['significant_dataset_count'])}",
            va="center",
            fontsize=7,
        )
    fig.tight_layout()
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(MECH_DIR / f"dhx_ddx_helicase_priority.{suffix}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_string_network_plot(edges: pd.DataFrame, nodes: pd.DataFrame) -> None:
    if edges.empty or nodes.empty or "preferredName_A" not in edges.columns:
        return
    node_names = nodes["gene_symbol"].tolist()
    angle = np.linspace(0, 2 * np.pi, len(node_names), endpoint=False)
    positions = {name: (np.cos(a), np.sin(a)) for name, a in zip(node_names, angle, strict=True)}
    score_norm = colors.Normalize(vmin=-2.5, vmax=2.5)
    cmap = plt.get_cmap("coolwarm")
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_axis_off()
    for _, edge in edges.iterrows():
        a = str(edge["preferredName_A"])
        b = str(edge["preferredName_B"])
        if a not in positions or b not in positions:
            continue
        score = float(edge.get("score", edge.get("combined_score", 0.7)))
        ax.plot(
            [positions[a][0], positions[b][0]],
            [positions[a][1], positions[b][1]],
            color="#555555",
            alpha=max(0.12, min(score, 1.0) * 0.45),
            lw=0.5 + min(score, 1.0) * 2.0,
            zorder=1,
        )
    for _, node in nodes.iterrows():
        name = node["gene_symbol"]
        x, y = positions[name]
        lfc = pd.to_numeric(node.get("log2FoldChange"), errors="coerce")
        color = cmap(score_norm(lfc)) if pd.notna(lfc) else "#CCCCCC"
        degree = pd.to_numeric(node.get("degree_in_string_subnetwork"), errors="coerce")
        size = 120 + 45 * (float(degree) if pd.notna(degree) else 0.0)
        ax.scatter(x, y, s=size, color=color, edgecolor="#222222", linewidth=0.6, zorder=2)
        ax.text(x * 1.12, y * 1.12, name, ha="center", va="center", fontsize=8)
    ax.set_title("High-confidence STRING physical subnetwork for prioritized targets")
    fig.tight_layout()
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(MECH_DIR / f"string_priority_network.{suffix}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    MECH_DIR.mkdir(parents=True, exist_ok=True)
    annotated = load_annotated_results(args.padj, args.effect_size, args.force_gencode)
    annotated.to_csv(MECH_DIR / "dge_results_annotated_host_biotype.csv", index=False)

    curated = annotated[annotated["translation_target"]].copy()
    curated.to_csv(MECH_DIR / "translation_targets_biotype_audit.csv", index=False)

    load = viral_load_table()
    load.to_csv(MECH_DIR / "viral_load_by_sample.csv", index=False)

    correlations = correlate_host_targets_with_viral_load(annotated)
    correlations.to_csv(MECH_DIR / "host_gene_viral_load_correlations.csv", index=False)

    enrichment = module_enrichment(annotated, args.padj, args.effect_size)
    enrichment.to_csv(MECH_DIR / "curated_module_enrichment.csv", index=False)

    priority = build_priority_table(annotated, correlations, args.padj, args.effect_size).reset_index(drop=True)
    priority.insert(0, "priority_rank", np.arange(1, len(priority) + 1))
    priority.to_csv(MECH_DIR / "mechanistic_priority_rank.csv", index=False)
    priority[priority["target_family"].isin(["DHX", "DDX"])].to_csv(MECH_DIR / "dhx_ddx_helicase_priority_rank.csv", index=False)
    priority.head(args.top_n).to_csv(MECH_DIR / "top_mechanistic_candidates.csv", index=False)
    write_summary_metrics(annotated, priority, enrichment)

    string_symbols = priority.loc[
        priority["protein_coding"] & ~priority["symbol_pseudogene_like"], "gene_symbol"
    ].head(min(args.top_n, 30)).dropna().astype(str).tolist()
    edges = fetch_string_network(string_symbols, force=args.force_string, skip=args.skip_string)
    if not edges.empty:
        edges.to_csv(MECH_DIR / "string_network_edges.csv", index=False)
    nodes = build_string_nodes(edges, priority)

    save_host_curated_volcano(annotated, priority, args.padj, args.top_n)
    save_dhx29_viral_load_plot(priority)
    save_module_enrichment_plot(enrichment)
    save_candidate_evidence_matrix(priority, args.top_n)
    save_helicase_priority_plot(priority)
    save_string_network_plot(edges, nodes)

    dhx29 = priority[priority["gene_symbol"].eq("DHX29")]
    dhx29_text = "absent"
    if not dhx29.empty:
        row = dhx29.iloc[0]
        dhx29_text = f"rank={int(row['priority_rank'])}, score={row['candidate_priority_score']:.3f}"
    print(f"Curated host protein-coding translation targets: {len(priority)}")
    print(f"Top candidate: {priority.iloc[0]['gene_symbol']} score={priority.iloc[0]['candidate_priority_score']:.3f}")
    print(f"DHX29 mechanistic priority: {dhx29_text}")
    print(f"Wrote mechanistic evidence outputs to {MECH_DIR}")


if __name__ == "__main__":
    main()
