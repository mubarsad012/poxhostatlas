#!/usr/bin/env python3
"""Filter and align raw count tables for PyDESeq2."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
PROCESSED_DIR = REPO_ROOT / "data" / "processed"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--min-total-count",
        type=int,
        default=10,
        help="Minimum total reads across selected samples required to keep a gene.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    counts_path = PROCESSED_DIR / "counts.csv"
    metadata_path = PROCESSED_DIR / "metadata.csv"

    counts = pd.read_csv(counts_path, dtype={"gene_id": "string", "gene_symbol": "string"})
    metadata = pd.read_csv(metadata_path, index_col="sample_id")

    sample_columns = metadata.index.tolist()
    missing = [sample for sample in sample_columns if sample not in counts.columns]
    if missing:
        raise SystemExit(f"Counts table is missing metadata samples: {missing}")

    count_values = counts[sample_columns].apply(pd.to_numeric, downcast="integer")
    if (count_values < 0).any().any():
        raise SystemExit("Counts table contains negative values.")

    total_counts = count_values.sum(axis=1)
    keep_mask = total_counts >= args.min_total_count
    filtered_counts = counts.loc[keep_mask].copy()
    filtered_values = count_values.loc[keep_mask].astype("int64")

    gene_annotation = filtered_counts[["gene_id", "gene_symbol"]].copy()
    gene_annotation["total_count"] = total_counts.loc[keep_mask].to_numpy()

    model_ready_counts = filtered_values.copy()
    model_ready_counts.insert(0, "sample_id_placeholder", filtered_counts["gene_id"].to_numpy())
    model_ready_counts = model_ready_counts.set_index("sample_id_placeholder").T
    model_ready_counts.index.name = "sample_id"
    model_ready_counts.columns.name = None
    model_ready_counts = model_ready_counts.loc[sample_columns]
    model_ready_counts = model_ready_counts.astype("int64")

    if model_ready_counts.index.tolist() != metadata.index.tolist():
        raise SystemExit("Model-ready count rows do not align exactly with metadata rows.")

    model_metadata = metadata.copy()
    model_metadata["infection"] = pd.Categorical(
        model_metadata["infection"],
        categories=["mock", "VacV"],
        ordered=True,
    )

    model_ready_counts.to_csv(PROCESSED_DIR / "model_ready_counts.csv")
    model_metadata.to_csv(PROCESSED_DIR / "model_ready_metadata.csv")
    gene_annotation.to_csv(PROCESSED_DIR / "gene_annotation.csv", index=False)

    dropped = len(counts) - len(filtered_counts)
    print(f"Kept {len(filtered_counts)} genes; dropped {dropped} genes with total count < {args.min_total_count}.")
    print(f"Wrote {PROCESSED_DIR / 'model_ready_counts.csv'} and aligned metadata.")


if __name__ == "__main__":
    main()
