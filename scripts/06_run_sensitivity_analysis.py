#!/usr/bin/env python3
"""This is to run strict-title sensitivity analysis that is excluding metadata-conflict samples."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats


REPO_ROOT = Path(__file__).resolve().parents[1]
PROCESSED_DIR = REPO_ROOT / "data" / "processed"
TABLE_DIR = REPO_ROOT / "results" / "tables"
SENSITIVITY_DIR = REPO_ROOT / "results" / "sensitivity"
TARGET_PATTERN = re.compile(r"^(DHX|EIF|RPS|RPL)", re.IGNORECASE)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--min-total-count", type=int, default=10)
    parser.add_argument("--padj", type=float, default=0.05)
    parser.add_argument("--n-cpus", type=int, default=8)
    return parser.parse_args()


def prepare_strict_inputs(min_total_count: int) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    counts = pd.read_csv(PROCESSED_DIR / "counts.csv", dtype={"gene_id": "string", "gene_symbol": "string"})
    metadata = pd.read_csv(PROCESSED_DIR / "metadata.csv", index_col="sample_id")
    manifest = pd.read_csv(PROCESSED_DIR / "sample_manifest.csv")

    conflict_samples = manifest.loc[
        manifest["metadata_warning"].notna() & (manifest["metadata_warning"] != ""),
        "sample_id",
    ].tolist()
    strict_metadata = metadata.drop(index=conflict_samples)
    sample_columns = strict_metadata.index.tolist()

    strict_values = counts[sample_columns].apply(pd.to_numeric, downcast="integer")
    keep_mask = strict_values.sum(axis=1) >= min_total_count
    strict_counts = strict_values.loc[keep_mask].astype("int64").T
    strict_counts.index.name = "sample_id"
    strict_counts.columns = counts.loc[keep_mask, "gene_id"].to_numpy()
    strict_counts.columns.name = None
    strict_counts = strict_counts.loc[sample_columns]

    gene_annotation = counts.loc[keep_mask, ["gene_id", "gene_symbol"]].copy()
    gene_annotation["strict_total_count"] = strict_values.loc[keep_mask].sum(axis=1).to_numpy()

    return strict_counts, strict_metadata, gene_annotation


def run_deseq2(counts: pd.DataFrame, metadata: pd.DataFrame, gene_annotation: pd.DataFrame, n_cpus: int) -> pd.DataFrame:
    if counts.index.tolist() != metadata.index.tolist():
        raise SystemExit("Strict-title count rows and metadata rows do not align.")
    if metadata["infection"].value_counts().to_dict() != {"VacV": 5, "mock": 2}:
        raise SystemExit(f"Unexpected strict-title sample counts: {metadata['infection'].value_counts().to_dict()}")

    inference = DefaultInference(n_cpus=n_cpus)
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
    ordered = ["gene_id", "gene_symbol", *[col for col in results.columns if col not in {"gene_id", "gene_symbol"}]]
    return results[ordered]


def extract_translation(results: pd.DataFrame, padj: float) -> tuple[pd.DataFrame, pd.DataFrame]:
    matched = results[results["gene_symbol"].fillna("").str.match(TARGET_PATTERN)].copy()
    matched = matched.sort_values(["log2FoldChange", "padj"], ascending=[True, True])
    matched["regulation_direction"] = matched["log2FoldChange"].map(
        lambda value: "upregulated" if value > 0 else "downregulated" if value < 0 else "unchanged"
    )
    impact = matched[matched["padj"].notna() & (matched["padj"] < padj)].copy()
    impact = impact.sort_values("log2FoldChange")
    return matched, impact


def compare_primary_to_strict(strict_translation: pd.DataFrame) -> pd.DataFrame:
    primary = pd.read_csv(TABLE_DIR / "translation_factors_all.csv", dtype={"gene_id": "string", "gene_symbol": "string"})
    comparison = primary.merge(
        strict_translation,
        on=["gene_id", "gene_symbol"],
        how="inner",
        suffixes=("_primary", "_strict"),
    )
    comparison["direction_agreement"] = (
        comparison["log2FoldChange_primary"].fillna(0).map(lambda value: 1 if value > 0 else -1 if value < 0 else 0)
        == comparison["log2FoldChange_strict"].fillna(0).map(lambda value: 1 if value > 0 else -1 if value < 0 else 0)
    )
    comparison["significant_primary"] = comparison["padj_primary"].notna() & (comparison["padj_primary"] < 0.05)
    comparison["significant_strict"] = comparison["padj_strict"].notna() & (comparison["padj_strict"] < 0.05)
    return comparison.sort_values(["significant_primary", "significant_strict", "padj_primary"], ascending=[False, False, True])


def main() -> None:
    args = parse_args()
    SENSITIVITY_DIR.mkdir(parents=True, exist_ok=True)

    counts, metadata, gene_annotation = prepare_strict_inputs(args.min_total_count)
    counts.to_csv(SENSITIVITY_DIR / "model_ready_counts_strict_titles.csv")
    metadata.to_csv(SENSITIVITY_DIR / "model_ready_metadata_strict_titles.csv")
    gene_annotation.to_csv(SENSITIVITY_DIR / "gene_annotation_strict_titles.csv", index=False)

    results = run_deseq2(counts, metadata, gene_annotation, args.n_cpus)
    results.to_csv(SENSITIVITY_DIR / "dge_results_strict_titles.csv", index=False)

    translation_all, translation_impact = extract_translation(results, args.padj)
    translation_all.to_csv(SENSITIVITY_DIR / "translation_factors_all_strict_titles.csv", index=False)
    translation_impact.to_csv(SENSITIVITY_DIR / "translation_factors_impact_strict_titles.csv", index=False)

    comparison = compare_primary_to_strict(translation_all)
    comparison.to_csv(SENSITIVITY_DIR / "sensitivity_primary_vs_strict_translation.csv", index=False)

    dhx29 = comparison[comparison["gene_symbol"].eq("DHX29")]
    print(f"Strict-title samples: {len(metadata)}; infection counts: {metadata['infection'].value_counts().to_dict()}")
    print(f"Strict-title genes modeled: {counts.shape[1]}")
    print(f"Strict-title translation matches: {len(translation_all)}; significant at padj < {args.padj}: {len(translation_impact)}")
    if not dhx29.empty:
        row = dhx29.iloc[0]
        print(
            "DHX29 primary vs strict: "
            f"log2FC {row['log2FoldChange_primary']:.3f} -> {row['log2FoldChange_strict']:.3f}; "
            f"padj {row['padj_primary']:.3g} -> {row['padj_strict']:.3g}"
        )
    print(f"Wrote strict-title sensitivity outputs to {SENSITIVITY_DIR}.")


if __name__ == "__main__":
    main()
