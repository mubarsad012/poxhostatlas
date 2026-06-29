#!/usr/bin/env python3
"""Extract translation-initiation and ribosomal factors from DGE results."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
TABLE_DIR = REPO_ROOT / "results" / "tables"
TARGET_PATTERN = re.compile(r"^(DHX|EIF|RPS|RPL)", re.IGNORECASE)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--padj", type=float, default=0.05, help="Adjusted p-value threshold.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    results = pd.read_csv(TABLE_DIR / "dge_results_full.csv", dtype={"gene_id": "string", "gene_symbol": "string"})

    matched = results[results["gene_symbol"].fillna("").str.match(TARGET_PATTERN)].copy()
    matched = matched.sort_values(["log2FoldChange", "padj"], ascending=[True, True])
    matched["regulation_direction"] = matched["log2FoldChange"].map(
        lambda value: "upregulated" if value > 0 else "downregulated" if value < 0 else "unchanged"
    )

    significant = matched[matched["padj"].notna() & (matched["padj"] < args.padj)].copy()
    significant = significant.sort_values("log2FoldChange", ascending=True)

    all_path = TABLE_DIR / "translation_factors_all.csv"
    impact_path = TABLE_DIR / "translation_factors_impact.csv"
    matched.to_csv(all_path, index=False)
    significant.to_csv(impact_path, index=False)

    print(f"Matched {len(matched)} DHX/EIF/RPS/RPL genes; {len(significant)} meet padj < {args.padj}.")
    print(f"Wrote {all_path} and {impact_path}.")


if __name__ == "__main__":
    main()
