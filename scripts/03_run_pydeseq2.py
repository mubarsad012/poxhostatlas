#!/usr/bin/env python3
"""Run PyDESeq2 differential expression for VacV versus mock."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats


REPO_ROOT = Path(__file__).resolve().parents[1]
PROCESSED_DIR = REPO_ROOT / "data" / "processed"
TABLE_DIR = REPO_ROOT / "results" / "tables"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-cpus", type=int, default=8, help="Number of CPUs for PyDESeq2 inference.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    TABLE_DIR.mkdir(parents=True, exist_ok=True)

    counts = pd.read_csv(PROCESSED_DIR / "model_ready_counts.csv", index_col="sample_id")
    metadata = pd.read_csv(PROCESSED_DIR / "model_ready_metadata.csv", index_col="sample_id")
    gene_annotation = pd.read_csv(PROCESSED_DIR / "gene_annotation.csv", dtype={"gene_id": "string", "gene_symbol": "string"})

    if counts.index.tolist() != metadata.index.tolist():
        raise SystemExit("Count rows and metadata rows must align exactly before PyDESeq2.")
    if not counts.apply(lambda col: pd.api.types.is_integer_dtype(col)).all():
        counts = counts.astype("int64")

    inference = DefaultInference(n_cpus=args.n_cpus)
    dds = DeseqDataSet(
        counts=counts,
        metadata=metadata,
        design="~infection",
        refit_cooks=True,
        inference=inference,
    )
    dds.deseq2()

    stats = DeseqStats(dds, contrast=["infection", "VacV", "mock"], inference=inference)
    stats.summary()

    results = stats.results_df.copy()
    results.index.name = "gene_id"
    results = results.reset_index()
    results = results.merge(gene_annotation[["gene_id", "gene_symbol"]], on="gene_id", how="left")
    ordered_columns = ["gene_id", "gene_symbol", *[col for col in results.columns if col not in {"gene_id", "gene_symbol"}]]
    results = results[ordered_columns]
    results.to_csv(TABLE_DIR / "dge_results_full.csv", index=False)

    normalized_counts = pd.DataFrame(
        dds.layers["normed_counts"],
        index=counts.index,
        columns=counts.columns,
    )
    normalized_counts.to_csv(PROCESSED_DIR / "normalized_counts.csv")

    required = {"log2FoldChange", "pvalue", "padj"}
    missing = required - set(results.columns)
    if missing:
        raise SystemExit(f"PyDESeq2 results missing expected columns: {sorted(missing)}")

    print(f"Wrote {TABLE_DIR / 'dge_results_full.csv'} with {len(results)} genes.")
    print(f"Wrote {PROCESSED_DIR / 'normalized_counts.csv'}.")


if __name__ == "__main__":
    main()
