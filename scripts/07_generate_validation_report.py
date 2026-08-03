#!/usr/bin/env python3
"""Generate validation report, key result tables, and hash receipt."""

from __future__ import annotations

import hashlib
import platform
import subprocess
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
EXTERNAL_RAW_DIR = REPO_ROOT / "data" / "external" / "expansion_raw"
PROCESSED_DIR = REPO_ROOT / "data" / "processed"
RESULTS_DIR = REPO_ROOT / "results"
TABLE_DIR = RESULTS_DIR / "tables"
FIGURE_DIR = RESULTS_DIR / "figures"
SENSITIVITY_DIR = RESULTS_DIR / "sensitivity"
META_DIR = RESULTS_DIR / "meta"
MECH_DIR = RESULTS_DIR / "mechanistic"
REPORT_DIR = RESULTS_DIR / "reports"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_value(*args: str) -> str:
    try:
        return subprocess.check_output(["git", *args], cwd=REPO_ROOT, text=True).strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "unavailable"


def metric_value(metrics: dict[str, object], name: str) -> str:
    value = metrics.get(name, "unavailable")
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return str(value)
    if numeric.is_integer():
        return str(int(numeric))
    return f"{numeric:.6g}"


def tracked_and_generated_files() -> list[Path]:
    candidates = [
        REPO_ROOT / "README.md",
        REPO_ROOT / "requirements.txt",
        REPO_ROOT / "requirements-lock.txt",
        REPO_ROOT / "docs" / "NOVELTY_AUDIT.md",
        REPO_ROOT / "docs" / "EXPANSION_DATASET_REGISTRY.csv",
        REPO_ROOT / "docs" / "public_dataset_discovery_results.csv",
        PROCESSED_DIR / "sample_manifest.csv",
        PROCESSED_DIR / "counts.csv",
        PROCESSED_DIR / "metadata.csv",
        PROCESSED_DIR / "model_ready_counts.csv",
        PROCESSED_DIR / "model_ready_metadata.csv",
        PROCESSED_DIR / "normalized_counts.csv",
        PROCESSED_DIR / "gencode_v44_gene_biotypes.csv",
        TABLE_DIR / "dge_results_full.csv",
        TABLE_DIR / "translation_factors_all.csv",
        TABLE_DIR / "translation_factors_impact.csv",
        SENSITIVITY_DIR / "dge_results_strict_titles.csv",
        SENSITIVITY_DIR / "translation_factors_all_strict_titles.csv",
        SENSITIVITY_DIR / "translation_factors_impact_strict_titles.csv",
        SENSITIVITY_DIR / "sensitivity_primary_vs_strict_translation.csv",
        FIGURE_DIR / "volcano_plot.png",
        FIGURE_DIR / "volcano_plot.pdf",
        FIGURE_DIR / "volcano_plot.svg",
        FIGURE_DIR / "translation_heatmap.png",
        FIGURE_DIR / "translation_heatmap.pdf",
        FIGURE_DIR / "translation_heatmap.svg",
        META_DIR / "cross_dataset_pilot_registry.csv",
        META_DIR / "GSE287860_translation_targets.csv",
        META_DIR / "GSE288000_NTC_translation_targets.csv",
        META_DIR / "poxvirus_translation_factor_meta_signature.csv",
        META_DIR / "poxvirus_translation_factor_meta_signature.png",
        META_DIR / "poxvirus_translation_factor_meta_signature.pdf",
        META_DIR / "poxvirus_translation_factor_meta_signature.svg",
        MECH_DIR / "mechanistic_summary_metrics.csv",
        MECH_DIR / "mechanistic_priority_rank.csv",
        MECH_DIR / "dhx_ddx_helicase_priority_rank.csv",
        MECH_DIR / "top_mechanistic_candidates.csv",
        MECH_DIR / "curated_module_enrichment.csv",
        MECH_DIR / "host_gene_viral_load_correlations.csv",
        MECH_DIR / "viral_load_by_sample.csv",
        MECH_DIR / "translation_targets_biotype_audit.csv",
        MECH_DIR / "string_network_edges.csv",
        MECH_DIR / "string_network_nodes.csv",
        MECH_DIR / "host_only_curated_volcano.png",
        MECH_DIR / "host_only_curated_volcano.pdf",
        MECH_DIR / "host_only_curated_volcano.svg",
        MECH_DIR / "candidate_evidence_matrix.png",
        MECH_DIR / "candidate_evidence_matrix.pdf",
        MECH_DIR / "candidate_evidence_matrix.svg",
        MECH_DIR / "dhx_ddx_helicase_priority.png",
        MECH_DIR / "dhx_ddx_helicase_priority.pdf",
        MECH_DIR / "dhx_ddx_helicase_priority.svg",
        MECH_DIR / "dhx29_viral_load_correlation.png",
        MECH_DIR / "dhx29_viral_load_correlation.pdf",
        MECH_DIR / "dhx29_viral_load_correlation.svg",
        MECH_DIR / "module_enrichment_barplot.png",
        MECH_DIR / "module_enrichment_barplot.pdf",
        MECH_DIR / "module_enrichment_barplot.svg",
        MECH_DIR / "string_priority_network.png",
        MECH_DIR / "string_priority_network.pdf",
        MECH_DIR / "string_priority_network.svg",
    ]
    candidates.extend(sorted((REPO_ROOT / "scripts").glob("*.py")))
    candidates.extend(sorted(EXTERNAL_RAW_DIR.glob("*.gz")))
    return [path for path in candidates if path.exists()]


def main() -> None:
    REPORT_DIR.mkdir(parents=True, exist_ok=True)

    metadata = pd.read_csv(PROCESSED_DIR / "model_ready_metadata.csv", index_col="sample_id")
    counts = pd.read_csv(PROCESSED_DIR / "model_ready_counts.csv", index_col="sample_id")
    manifest = pd.read_csv(PROCESSED_DIR / "sample_manifest.csv")
    results = pd.read_csv(TABLE_DIR / "dge_results_full.csv")
    translation_all = pd.read_csv(TABLE_DIR / "translation_factors_all.csv")
    translation_impact = pd.read_csv(TABLE_DIR / "translation_factors_impact.csv")
    strict_metadata = pd.read_csv(SENSITIVITY_DIR / "model_ready_metadata_strict_titles.csv", index_col="sample_id")
    strict_counts = pd.read_csv(SENSITIVITY_DIR / "model_ready_counts_strict_titles.csv", index_col="sample_id")
    strict_impact = pd.read_csv(SENSITIVITY_DIR / "translation_factors_impact_strict_titles.csv")
    comparison = pd.read_csv(SENSITIVITY_DIR / "sensitivity_primary_vs_strict_translation.csv")
    meta_signature_path = META_DIR / "poxvirus_translation_factor_meta_signature.csv"
    meta_registry_path = META_DIR / "cross_dataset_pilot_registry.csv"
    meta_signature = pd.read_csv(meta_signature_path) if meta_signature_path.exists() else pd.DataFrame()
    meta_registry = pd.read_csv(meta_registry_path) if meta_registry_path.exists() else pd.DataFrame()
    mech_summary_path = MECH_DIR / "mechanistic_summary_metrics.csv"
    mech_priority_path = MECH_DIR / "mechanistic_priority_rank.csv"
    mech_enrichment_path = MECH_DIR / "curated_module_enrichment.csv"
    mech_summary = pd.read_csv(mech_summary_path) if mech_summary_path.exists() else pd.DataFrame()
    mech_priority = pd.read_csv(mech_priority_path) if mech_priority_path.exists() else pd.DataFrame()
    mech_enrichment = pd.read_csv(mech_enrichment_path) if mech_enrichment_path.exists() else pd.DataFrame()

    dhx29 = results[results["gene_symbol"].eq("DHX29")].iloc[0]
    strict_dhx29 = comparison[comparison["gene_symbol"].eq("DHX29")].iloc[0]
    metadata_warnings = manifest[
        manifest["metadata_warning"].notna() & (manifest["metadata_warning"] != "")
    ][["sample_id", "metadata_warning"]]

    top_down = translation_impact.sort_values("log2FoldChange").head(10)
    top_up = translation_impact.sort_values("log2FoldChange", ascending=False).head(10)
    top_down.to_csv(TABLE_DIR / "top_downregulated_translation_factors.csv", index=False)
    top_up.to_csv(TABLE_DIR / "top_upregulated_translation_factors.csv", index=False)

    hash_rows = []
    for path in tracked_and_generated_files():
        hash_rows.append(
            {
                "relative_path": path.relative_to(REPO_ROOT).as_posix(),
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    hash_table = pd.DataFrame(hash_rows).sort_values("relative_path")
    hash_table.to_csv(REPORT_DIR / "file_hashes_sha256.csv", index=False)

    direction_subset = comparison[comparison["significant_primary"] | comparison["significant_strict"]]
    direction_rate = direction_subset["direction_agreement"].mean() if not direction_subset.empty else float("nan")

    lines = [
        "# Validation Report",
        "",
        "## Execution Context",
        "",
        f"- Git commit: `{git_value('rev-parse', '--short', 'HEAD')}`",
        f"- Git status: `{git_value('status', '--short') or 'clean except ignored local artifacts'}`",
        f"- Python: `{platform.python_version()}`",
        f"- Platform: `{platform.platform()}`",
        "",
        "## Primary Model",
        "",
        f"- Samples: {len(metadata)} total.",
        f"- Infection counts: {metadata['infection'].value_counts().to_dict()}.",
        f"- Model-ready genes: {counts.shape[1]}.",
        f"- Count/metadata alignment: {counts.index.tolist() == metadata.index.tolist()}.",
        f"- Full DGE rows: {len(results)}.",
        f"- Required DGE columns present: { {'log2FoldChange', 'pvalue', 'padj'}.issubset(results.columns) }.",
        f"- Translation-factor matches: {len(translation_all)}.",
        f"- Significant translation-factor impacts at padj < 0.05: {len(translation_impact)}.",
        "",
        "## DHX29 Primary Signal",
        "",
        f"- Gene ID: `{dhx29['gene_id']}`.",
        f"- log2 fold change, VacV versus mock: {dhx29['log2FoldChange']:.6f}.",
        f"- Adjusted p-value: {dhx29['padj']:.6g}.",
        "",
        "## Metadata Conflict Handling",
        "",
        "The primary model uses supplementary count filenames as the fraction source of truth. The following GEO title/file conflicts were recorded and retained in the manifest:",
        "",
    ]
    for _, row in metadata_warnings.iterrows():
        lines.append(f"- `{row['sample_id']}`: `{row['metadata_warning']}`")
    lines.extend(
        [
            "",
            "## Strict-Title Sensitivity Model",
            "",
            f"- Excluded samples: {metadata_warnings['sample_id'].tolist()}.",
            f"- Samples: {len(strict_metadata)} total.",
            f"- Infection counts: {strict_metadata['infection'].value_counts().to_dict()}.",
            f"- Model-ready genes: {strict_counts.shape[1]}.",
            f"- Significant strict-title translation-factor impacts at padj < 0.05: {len(strict_impact)}.",
            f"- DHX29 primary log2FC/padj: {strict_dhx29['log2FoldChange_primary']:.6f} / {strict_dhx29['padj_primary']:.6g}.",
            f"- DHX29 strict-title log2FC/padj: {strict_dhx29['log2FoldChange_strict']:.6f} / {strict_dhx29['padj_strict']:.6g}.",
            f"- Direction agreement among translation factors significant in either primary or strict model: {direction_rate:.3f}.",
            "",
            "## Cross-Dataset Novelty Pilot",
            "",
        ]
    )
    if meta_signature.empty:
        lines.append("- Cross-dataset pilot outputs were not present when this report was generated.")
    else:
        lines.extend(
            [
                f"- Meta-signature rows: {len(meta_signature)}.",
                f"- Pilot registry rows: {len(meta_registry)}.",
            ]
        )
        for _, row in meta_registry.iterrows():
            lines.append(
                "- `{dataset}`: {genes_modeled} genes modeled; {target_genes} target genes; "
                "{significant_targets} significant targets.".format(**row.to_dict())
            )
        dhx29_meta = meta_signature[meta_signature["gene_symbol"].eq("DHX29")]
        if not dhx29_meta.empty:
            row = dhx29_meta.iloc[0]
            lines.extend(
                [
                    f"- DHX29 dataset count: {int(row['dataset_count'])}.",
                    f"- DHX29 significant dataset count: {int(row['significant_dataset_count'])}.",
                    f"- DHX29 mean log2 fold change across pilot contexts: {row['mean_log2_fold_change']:.6f}.",
                    f"- DHX29 direction consistency: {row['direction_consistency']:.3f}.",
                ]
            )
    lines.extend(["", "## Mechanistic Evidence Layer", ""])
    if mech_summary.empty:
        lines.append("- Mechanistic prioritization outputs were not present when this report was generated.")
    else:
        metrics = dict(zip(mech_summary["metric"], mech_summary["value"], strict=False))
        lines.extend(
            [
                f"- Host model genes: {metric_value(metrics, 'host_model_genes')}.",
                f"- Viral model genes excluded from host-only curation: {metric_value(metrics, 'viral_model_genes')}.",
                f"- Host protein-coding model genes by GENCODE v44: {metric_value(metrics, 'host_protein_coding_model_genes')}.",
                f"- Curated host protein-coding translation targets: {metric_value(metrics, 'curated_host_protein_coding_translation_targets')}.",
                f"- Curated significant targets with abs(log2FC) >= 0.5: {metric_value(metrics, 'curated_significant_with_abs_log2fc_ge_threshold')}.",
                f"- DHX29 overall mechanistic priority rank: {metric_value(metrics, 'dhx29_priority_rank')}.",
                f"- DHX29 DHX/DDX helicase priority rank: {metric_value(metrics, 'dhx29_helicase_priority_rank')}.",
                f"- Minimum curated module enrichment padj: {metric_value(metrics, 'module_enrichment_min_padj')}.",
            ]
        )
        if not mech_priority.empty:
            top = mech_priority.iloc[0]
            lines.append(
                f"- Top mechanistic candidate: `{top['gene_symbol']}` "
                f"(score {top['candidate_priority_score']:.3f}; log2FC {top['log2FoldChange']:.3f}; padj {top['padj']:.3g})."
            )
            dhx29_priority = mech_priority[mech_priority["gene_symbol"].eq("DHX29")]
            if not dhx29_priority.empty:
                row = dhx29_priority.iloc[0]
                lines.append(
                    f"- DHX29 viral-load correlation in infected samples: Pearson r={row['pearson_r_viral_load']:.3f}, "
                    f"p={row['pearson_p_viral_load']:.3g}; exploratory because n={int(row['n_infected_samples'])} infected samples."
                )
        if not mech_enrichment.empty:
            enrich = mech_enrichment.iloc[0]
            lines.append(
                f"- Strongest curated module enrichment: `{enrich['query_set']}` / `{enrich['module']}` "
                f"(overlap {int(enrich['overlap_size'])}; padj {enrich['padj']:.3g})."
            )
    lines.extend(["", "## Figure Checks", ""])
    for figure in sorted(FIGURE_DIR.glob("*")):
        if figure.suffix.lower() in {".png", ".pdf", ".svg"}:
            lines.append(f"- `{figure.relative_to(REPO_ROOT)}`: {figure.stat().st_size} bytes; SHA256 `{sha256_file(figure)}`.")
    for figure in sorted(META_DIR.glob("*")):
        if figure.suffix.lower() in {".png", ".pdf", ".svg"}:
            lines.append(f"- `{figure.relative_to(REPO_ROOT)}`: {figure.stat().st_size} bytes; SHA256 `{sha256_file(figure)}`.")
    for figure in sorted(MECH_DIR.glob("*")):
        if figure.suffix.lower() in {".png", ".pdf", ".svg"}:
            lines.append(f"- `{figure.relative_to(REPO_ROOT)}`: {figure.stat().st_size} bytes; SHA256 `{sha256_file(figure)}`.")
    lines.extend(
        [
            "",
            "## Hash Receipt",
            "",
            "A machine-readable hash table was written to `results/reports/file_hashes_sha256.csv`.",
            "",
            "## Interpretation Boundary",
            "",
            "These outputs support a transcript-level secondary analysis of a public dataset. They do not establish direct protein-level mechanism without additional translational, proteomic, or perturbation validation.",
            "",
        ]
    )

    (REPORT_DIR / "validation_report.md").write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote {REPORT_DIR / 'validation_report.md'}")
    print(f"Wrote {REPORT_DIR / 'file_hashes_sha256.csv'}")
    print(f"DHX29 primary log2FC={dhx29['log2FoldChange']:.3f}, padj={dhx29['padj']:.3g}")


if __name__ == "__main__":
    main()
