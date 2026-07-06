#!/usr/bin/env python3
"""Run pilot cross-dataset poxvirus translation-factor meta-signature."""

from __future__ import annotations

import argparse
import gzip
import re
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import requests
import seaborn as sns
from matplotlib import pyplot as plt
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats


REPO_ROOT = Path(__file__).resolve().parents[1]
EXPANSION_RAW = REPO_ROOT / "data" / "external" / "expansion_raw"
META_DIR = REPO_ROOT / "results" / "meta"
PROCESSED_DIR = REPO_ROOT / "data" / "processed"
TABLE_DIR = REPO_ROOT / "results" / "tables"
TARGET_PATTERN = re.compile(r"^(DHX|DDX|EIF|RPS|RPL)", re.IGNORECASE)


DATASETS = {
    "GSE287860": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE287nnn/GSE287860/suppl/GSE287860_raw_feature_counts.txt.gz",
        "filename": "GSE287860_raw_feature_counts.txt.gz",
        "columns": {
            "cherry_rep1": "control",
            "cherry_rep2": "control",
            "cherry_rep3": "control",
            "m3_rep1": "M003",
            "m3_rep2": "M003",
            "m3_rep3": "M003",
        },
        "contrast": ["condition", "M003", "control"],
        "description": "M003.1 versus mCherry compact poxvirus effector RNA-seq contrast",
    },
    "GSE288000_NTC": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE288nnn/GSE288000/suppl/GSE288000_raw_feature_counts2.txt.gz",
        "filename": "GSE288000_raw_feature_counts2.txt.gz",
        "columns": {
            "NTC_cherry_rep1": "control",
            "NTC_cherry_rep2": "control",
            "NTC_cherry_rep3": "control",
            "NTC_M003_rep1": "M003",
            "NTC_M003_rep2": "M003",
            "NTC_M003_rep3": "M003",
        },
        "contrast": ["condition", "M003", "control"],
        "description": "NTC background M003.1 versus mCherry poxvirus effector RNA-seq contrast",
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--force", action="store_true", help="Redownload expansion matrices.")
    parser.add_argument("--min-total-count", type=int, default=10)
    parser.add_argument("--padj", type=float, default=0.05)
    parser.add_argument("--n-cpus", type=int, default=8)
    return parser.parse_args()


def download(url: str, destination: Path, force: bool = False) -> None:
    if destination.exists() and destination.stat().st_size > 0 and not force:
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    tmp = destination.with_suffix(destination.suffix + ".tmp")
    with requests.get(url, stream=True, timeout=120) as response:
        response.raise_for_status()
        with tmp.open("wb") as handle:
            for chunk in response.iter_content(1024 * 1024):
                if chunk:
                    handle.write(chunk)
    tmp.replace(destination)


def gene_symbol_map() -> dict[str, str]:
    counts = pd.read_csv(PROCESSED_DIR / "counts.csv", usecols=["gene_id", "gene_symbol"], dtype="string")
    mapping: dict[str, str] = {}
    for _, row in counts.iterrows():
        gene_id = str(row["gene_id"])
        symbol = str(row["gene_symbol"])
        mapping[gene_id] = symbol
        mapping[gene_id.split(".", 1)[0]] = symbol
    return mapping


def read_count_matrix(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        frame = pd.read_csv(handle, sep="\t", index_col=0)
    frame.index = frame.index.astype(str).str.strip('"')
    frame.columns = frame.columns.astype(str).str.strip('"')
    return frame


def prepare_inputs(frame: pd.DataFrame, selected_columns: dict[str, str], min_total_count: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    missing = [column for column in selected_columns if column not in frame.columns]
    if missing:
        raise SystemExit(f"Missing expected expansion columns: {missing}")
    counts = frame[list(selected_columns)].apply(pd.to_numeric, errors="coerce").fillna(0).astype("int64")
    keep = counts.sum(axis=1) >= min_total_count
    counts = counts.loc[keep].T
    counts.index.name = "sample_id"
    metadata = pd.DataFrame(
        {
            "sample_id": list(selected_columns),
            "condition": [selected_columns[column] for column in selected_columns],
        }
    ).set_index("sample_id")
    return counts, metadata


def run_deseq2(counts: pd.DataFrame, metadata: pd.DataFrame, contrast: list[str], n_cpus: int) -> pd.DataFrame:
    inference = DefaultInference(n_cpus=n_cpus)
    dds = DeseqDataSet(
        counts=counts,
        metadata=metadata,
        design="~condition",
        refit_cooks=True,
        inference=inference,
    )
    dds.deseq2()
    stats = DeseqStats(dds, contrast=contrast, inference=inference)
    stats.summary()
    results = stats.results_df.copy()
    results.index.name = "gene_id"
    return results.reset_index()


def add_symbols_and_targets(results: pd.DataFrame, mapping: dict[str, str], padj: float) -> tuple[pd.DataFrame, pd.DataFrame]:
    results = results.copy()
    results["gene_id_base"] = results["gene_id"].astype(str).str.split(".", n=1).str[0]
    results["gene_symbol"] = results["gene_id"].map(mapping).fillna(results["gene_id_base"].map(mapping)).fillna(results["gene_id"])
    ordered = ["gene_id", "gene_symbol", *[col for col in results.columns if col not in {"gene_id", "gene_symbol"}]]
    results = results[ordered]
    targets = results[results["gene_symbol"].fillna("").str.match(TARGET_PATTERN)].copy()
    targets["significant"] = targets["padj"].notna() & (targets["padj"] < padj)
    targets["direction"] = np.where(targets["log2FoldChange"] > 0, "up", np.where(targets["log2FoldChange"] < 0, "down", "zero"))
    targets = targets.sort_values(["significant", "padj"], ascending=[False, True])
    return results, targets


def primary_target_frame() -> pd.DataFrame:
    primary = pd.read_csv(TABLE_DIR / "translation_factors_all.csv", dtype={"gene_id": "string", "gene_symbol": "string"})
    primary = primary[primary["gene_symbol"].fillna("").str.match(TARGET_PATTERN)].copy()
    primary["dataset"] = "GSE278320_primary"
    return primary[["dataset", "gene_id", "gene_symbol", "log2FoldChange", "padj"]]


def build_meta_signature(frames: list[pd.DataFrame], padj: float) -> pd.DataFrame:
    combined = pd.concat(frames, ignore_index=True)
    combined["significant"] = combined["padj"].notna() & (combined["padj"] < padj)
    pivot_lfc = combined.pivot_table(index="gene_symbol", columns="dataset", values="log2FoldChange", aggfunc="mean")
    pivot_padj = combined.pivot_table(index="gene_symbol", columns="dataset", values="padj", aggfunc="min")
    score = pivot_lfc.notna().sum(axis=1).rename("dataset_count").to_frame()
    score["significant_dataset_count"] = pivot_padj.lt(padj).sum(axis=1)
    score["mean_log2_fold_change"] = pivot_lfc.mean(axis=1)
    score["direction_consistency"] = pivot_lfc.apply(
        lambda row: max((row.dropna() > 0).mean(), (row.dropna() < 0).mean()) if row.notna().any() else np.nan,
        axis=1,
    )
    meta = score.join(pivot_lfc.add_prefix("log2fc_")).join(pivot_padj.add_prefix("padj_"))
    return meta.reset_index().sort_values(
        ["significant_dataset_count", "dataset_count", "direction_consistency"],
        ascending=[False, False, False],
    )


def write_meta_heatmap(meta: pd.DataFrame) -> None:
    lfc_cols = [col for col in meta.columns if col.startswith("log2fc_")]
    plot = meta[meta["significant_dataset_count"] > 0].head(40).set_index("gene_symbol")[lfc_cols]
    if plot.empty:
        return
    plot.columns = [col.replace("log2fc_", "") for col in plot.columns]
    sns.set_theme(style="white", context="paper")
    fig, ax = plt.subplots(figsize=(8, max(6, len(plot) * 0.18)))
    sns.heatmap(plot, cmap="vlag", center=0, linewidths=0.2, linecolor="#DDDDDD", ax=ax, cbar_kws={"label": "log2 fold change"})
    ax.set_title("Pilot cross-dataset poxvirus translation-factor signature")
    ax.set_xlabel("")
    ax.set_ylabel("")
    plt.tight_layout()
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(META_DIR / f"poxvirus_translation_factor_meta_signature.{suffix}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    EXPANSION_RAW.mkdir(parents=True, exist_ok=True)
    META_DIR.mkdir(parents=True, exist_ok=True)
    mapping = gene_symbol_map()

    meta_frames = [primary_target_frame()]
    registry_rows = []
    for dataset, spec in DATASETS.items():
        path = EXPANSION_RAW / spec["filename"]
        download(spec["url"], path, force=args.force)
        matrix = read_count_matrix(path)
        counts, metadata = prepare_inputs(matrix, spec["columns"], args.min_total_count)
        counts.to_csv(META_DIR / f"{dataset}_model_ready_counts.csv")
        metadata.to_csv(META_DIR / f"{dataset}_metadata.csv")
        results = run_deseq2(counts, metadata, spec["contrast"], args.n_cpus)
        results, targets = add_symbols_and_targets(results, mapping, args.padj)
        results.to_csv(META_DIR / f"{dataset}_dge_results.csv", index=False)
        targets.to_csv(META_DIR / f"{dataset}_translation_targets.csv", index=False)
        targets_for_meta = targets[["gene_id", "gene_symbol", "log2FoldChange", "padj"]].copy()
        targets_for_meta["dataset"] = dataset
        meta_frames.append(targets_for_meta[["dataset", "gene_id", "gene_symbol", "log2FoldChange", "padj"]])
        registry_rows.append(
            {
                "dataset": dataset,
                "description": spec["description"],
                "samples": len(metadata),
                "genes_modeled": counts.shape[1],
                "target_genes": len(targets),
                "significant_targets": int(targets["significant"].sum()),
                "source_url": spec["url"],
            }
        )
        print(f"{dataset}: {counts.shape[1]} genes modeled; {int(targets['significant'].sum())} significant targets.")

    meta = build_meta_signature(meta_frames, args.padj)
    meta.to_csv(META_DIR / "poxvirus_translation_factor_meta_signature.csv", index=False)
    pd.DataFrame(registry_rows).to_csv(META_DIR / "cross_dataset_pilot_registry.csv", index=False)
    write_meta_heatmap(meta)
    print(f"Wrote pilot meta-signature outputs to {META_DIR}.")


if __name__ == "__main__":
    main()
