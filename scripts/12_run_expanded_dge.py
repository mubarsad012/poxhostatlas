#!/usr/bin/env python3
"""Expand the poxvirus secondary analysis with additional independent contrasts.

This module multiplies the evidence base in two ways:

1. GSE288000 host-background panel: the GSE288000 feature-count matrix contains a
   five-way host genetic-background panel (NTC, N4BP1, ZC3H12A, TRIM25, and the
   ZC3H12A/N4BP1 double background). For each background we fit an independent
   PyDESeq2 model of the poxvirus M003 effector versus the mCherry control. This
   converts a single contrast into a controlled robustness panel that asks whether
   translation-factor remodeling is stable across antiviral host backgrounds.

2. GSE284044 Vaccinia time course: an independent Vaccinia-infected Vero data set
   with author-provided DESeq2 results at 2, 6, and 24 hours post infection. We
   harmonize these external result tables onto a common schema keyed by gene
   symbol so they can join the meta-analysis as a second Vaccinia study with
   temporal resolution.

All harmonized per-contrast tables are written to ``results/expanded`` and a
machine-readable registry of every contrast used downstream is emitted.
"""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

REPO_ROOT = Path(__file__).resolve().parents[1]
EXPANSION_RAW = REPO_ROOT / "data" / "external" / "expansion_raw"
PROCESSED_DIR = REPO_ROOT / "data" / "processed"
META_DIR = REPO_ROOT / "results" / "meta"
EXPANDED_DIR = REPO_ROOT / "results" / "expanded"

TARGET_PATTERN = re.compile(r"^(DHX|DDX|EIF|RPS|RPL)", re.IGNORECASE)
GSE288000_FILE = EXPANSION_RAW / "GSE288000_raw_feature_counts2.txt.gz"
GSE288000_BACKGROUNDS = ["NTC", "N4BP1", "ZC3H12A", "TRIM25", "ZC3H12A_N4BP1"]
GSE284044_TIMEPOINTS = [2, 6, 24]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--min-total-count", type=int, default=10)
    parser.add_argument("--padj", type=float, default=0.05)
    parser.add_argument("--n-cpus", type=int, default=4)
    return parser.parse_args()


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


def harmonize(results: pd.DataFrame, mapping: dict[str, str] | None, padj: float, symbol_col: str | None = None) -> pd.DataFrame:
    results = results.copy()
    if symbol_col is not None:
        results["gene_symbol"] = results[symbol_col].astype(str)
        results["gene_id"] = results[symbol_col].astype(str)
    else:
        results["gene_id_base"] = results["gene_id"].astype(str).str.split(".", n=1).str[0]
        results["gene_symbol"] = (
            results["gene_id"].map(mapping).fillna(results["gene_id_base"].map(mapping)).fillna(results["gene_id"])
        )
    keep = ["gene_id", "gene_symbol", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]
    keep = [c for c in keep if c in results.columns]
    results = results[keep]
    results["significant"] = results["padj"].notna() & (results["padj"] < padj)
    results["is_translation_factor"] = results["gene_symbol"].fillna("").str.match(TARGET_PATTERN)
    return results


def run_gse288000_panel(mapping: dict[str, str], min_total_count: int, padj: float, n_cpus: int) -> list[dict]:
    matrix = read_count_matrix(GSE288000_FILE)
    registry = []
    robustness_rows = []
    for background in GSE288000_BACKGROUNDS:
        cherry_cols = [f"{background}_cherry_rep{r}" for r in (1, 2, 3)]
        m003_cols = [f"{background}_M003_rep{r}" for r in (1, 2, 3)]
        selected = {c: "control" for c in cherry_cols}
        selected.update({c: "M003" for c in m003_cols})
        missing = [c for c in selected if c not in matrix.columns]
        if missing:
            raise SystemExit(f"GSE288000 background {background} missing columns: {missing}")
        counts = matrix[list(selected)].apply(pd.to_numeric, errors="coerce").fillna(0).astype("int64")
        counts = counts.loc[counts.sum(axis=1) >= min_total_count].T
        counts.index.name = "sample_id"
        metadata = pd.DataFrame(
            {"sample_id": list(selected), "condition": [selected[c] for c in selected]}
        ).set_index("sample_id")
        results = run_deseq2(counts, metadata, ["condition", "M003", "control"], n_cpus)
        results = harmonize(results, mapping, padj)
        contrast_id = f"GSE288000_{background}"
        results.to_csv(EXPANDED_DIR / f"{contrast_id}_dge.csv", index=False)
        tf = results[results["is_translation_factor"]]
        registry.append(
            {
                "contrast_id": contrast_id,
                "study": "GSE288000",
                "virus_system": "Myxoma_M003_effector",
                "context": f"{background} host background",
                "comparison": "M003 vs mCherry",
                "samples": int(counts.shape[0]),
                "genes_modeled": int(counts.shape[1]),
                "translation_factors": int(len(tf)),
                "significant_translation_factors": int(tf["significant"].sum()),
                "use_in_primary_meta": background == "NTC",
            }
        )
        for _, r in tf.iterrows():
            robustness_rows.append(
                {"gene_symbol": r["gene_symbol"], "background": background,
                 "log2FoldChange": r["log2FoldChange"], "padj": r["padj"]}
            )
        print(f"GSE288000 {background}: {counts.shape[1]} genes, {int(tf['significant'].sum())} sig TFs")
    robustness = pd.DataFrame(robustness_rows)
    if not robustness.empty:
        lfc = robustness.pivot_table(index="gene_symbol", columns="background", values="log2FoldChange", aggfunc="mean")
        lfc = lfc.reindex(columns=GSE288000_BACKGROUNDS)
        lfc["mean_log2FoldChange"] = lfc.mean(axis=1)
        lfc["direction_consistency"] = lfc[GSE288000_BACKGROUNDS].apply(
            lambda row: max((row.dropna() > 0).mean(), (row.dropna() < 0).mean()) if row.notna().any() else np.nan,
            axis=1,
        )
        lfc.sort_values("mean_log2FoldChange", ascending=False).to_csv(
            EXPANDED_DIR / "GSE288000_host_background_robustness.csv"
        )
    return registry


def run_gse284044_timecourse(padj: float) -> list[dict]:
    registry = []
    for hpi in GSE284044_TIMEPOINTS:
        path = EXPANSION_RAW / f"GSE284044_{hpi}hpi_DESeq2.csv.gz"
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            frame = pd.read_csv(handle, sep="\t")
        frame.columns = ["gene_symbol", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]
        results = harmonize(frame, None, padj, symbol_col="gene_symbol")
        contrast_id = f"GSE284044_{hpi}hpi"
        results.to_csv(EXPANDED_DIR / f"{contrast_id}_dge.csv", index=False)
        tf = results[results["is_translation_factor"]]
        registry.append(
            {
                "contrast_id": contrast_id,
                "study": "GSE284044",
                "virus_system": "Vaccinia_infection",
                "context": f"Vero {hpi} hpi",
                "comparison": "VacV vs mock",
                "samples": np.nan,
                "genes_modeled": int(len(results)),
                "translation_factors": int(len(tf)),
                "significant_translation_factors": int(tf["significant"].sum()),
                "use_in_primary_meta": hpi == 6,
            }
        )
        print(f"GSE284044 {hpi}hpi: {len(results)} genes, {int(tf['significant'].sum())} sig TFs")
    return registry


def main() -> None:
    args = parse_args()
    EXPANDED_DIR.mkdir(parents=True, exist_ok=True)
    mapping = gene_symbol_map()
    registry = []
    registry.extend(run_gse288000_panel(mapping, args.min_total_count, args.padj, args.n_cpus))
    registry.extend(run_gse284044_timecourse(args.padj))
    pd.DataFrame(registry).to_csv(EXPANDED_DIR / "expanded_contrast_registry.csv", index=False)
    print(f"\nWrote {len(registry)} expanded contrasts to {EXPANDED_DIR}")


if __name__ == "__main__":
    main()
