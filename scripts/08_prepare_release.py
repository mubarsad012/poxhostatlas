#!/usr/bin/env python3
"""Assemble manuscript-ready release package and archive."""

from __future__ import annotations

import hashlib
import shutil
import zipfile
from datetime import date
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
RELEASE_DIR = REPO_ROOT / "release"
PROCESSED_DIR = REPO_ROOT / "data" / "processed"
RAW_DIR = REPO_ROOT / "data" / "raw"
EXTERNAL_RAW_DIR = REPO_ROOT / "data" / "external" / "expansion_raw"
RESULTS_DIR = REPO_ROOT / "results"
TABLE_DIR = RESULTS_DIR / "tables"
FIGURE_DIR = RESULTS_DIR / "figures"
SENSITIVITY_DIR = RESULTS_DIR / "sensitivity"
META_DIR = RESULTS_DIR / "meta"
MECH_DIR = RESULTS_DIR / "mechanistic"
REPORT_DIR = RESULTS_DIR / "reports"
DOCS_DATASET_REGISTRY = REPO_ROOT / "docs" / "dataset_registry.csv"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def reset_release_tree() -> None:
    if RELEASE_DIR.exists():
        shutil.rmtree(RELEASE_DIR)
    for folder in [
        "manuscript",
        "figures",
        "tables",
        "data/raw/count_files",
        "data/raw/geo_metadata",
        "data/external/expansion_raw",
        "data/processed",
        "provenance",
        "code/scripts",
    ]:
        (RELEASE_DIR / folder).mkdir(parents=True, exist_ok=True)


def copy_file(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)


def write_pipeline_schematic() -> None:
    labels = [
        "GSE278320\nGEO metadata",
        "Official supplementary\ncount files",
        "Count matrix +\nmanifest QC",
        "PyDESeq2\nVacV vs mock",
        "Strict-title\nsensitivity",
        "GSE287860/GSE288000\npilot expansion",
        "GENCODE/STRING\nmechanistic layer",
        "Figures, tables,\nrelease archive",
    ]
    fig, ax = plt.subplots(figsize=(14, 3))
    ax.set_axis_off()
    x_positions = [0.055, 0.185, 0.315, 0.445, 0.575, 0.705, 0.835, 0.965]
    colors = ["#E8EEF7", "#EAF4EC", "#FFF2CC", "#FCE4D6", "#F4E4D6", "#EADCF8", "#DDEAF6", "#E2F0F3"]
    for idx, (x_pos, label, color) in enumerate(zip(x_positions, labels, colors, strict=True)):
        ax.text(
            x_pos,
            0.52,
            label,
            ha="center",
            va="center",
            fontsize=10,
            bbox={"boxstyle": "round,pad=0.35", "facecolor": color, "edgecolor": "#333333", "linewidth": 1.0},
        )
        if idx < len(x_positions) - 1:
            ax.annotate(
                "",
                xy=(x_positions[idx + 1] - 0.06, 0.52),
                xytext=(x_pos + 0.06, 0.52),
                arrowprops={"arrowstyle": "->", "lw": 1.2, "color": "#333333"},
            )
    fig.suptitle("Reproducible secondary-analysis workflow", fontsize=13, y=0.92)
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(RELEASE_DIR / "figures" / f"Figure_1_pipeline_overview.{suffix}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_manuscript() -> None:
    results = pd.read_csv(TABLE_DIR / "dge_results_full.csv")
    impact = pd.read_csv(TABLE_DIR / "translation_factors_impact.csv")
    comparison = pd.read_csv(SENSITIVITY_DIR / "sensitivity_primary_vs_strict_translation.csv")
    dhx29 = results[results["gene_symbol"].eq("DHX29")].iloc[0]
    strict_dhx29 = comparison[comparison["gene_symbol"].eq("DHX29")].iloc[0]
    top_eif3 = impact[impact["gene_symbol"].fillna("").str.startswith("EIF3")].sort_values("padj").head(5)
    meta_path = META_DIR / "poxvirus_translation_factor_meta_signature.csv"
    meta = pd.read_csv(meta_path) if meta_path.exists() else pd.DataFrame()
    dhx29_meta = meta[meta["gene_symbol"].eq("DHX29")].head(1) if not meta.empty else pd.DataFrame()
    mech_summary_path = MECH_DIR / "mechanistic_summary_metrics.csv"
    mech_priority_path = MECH_DIR / "mechanistic_priority_rank.csv"
    mech_enrichment_path = MECH_DIR / "curated_module_enrichment.csv"
    mech_summary = pd.read_csv(mech_summary_path) if mech_summary_path.exists() else pd.DataFrame()
    mech_priority = pd.read_csv(mech_priority_path) if mech_priority_path.exists() else pd.DataFrame()
    mech_enrichment = pd.read_csv(mech_enrichment_path) if mech_enrichment_path.exists() else pd.DataFrame()
    if dhx29_meta.empty:
        cross_dataset_text = "Cross-dataset pilot results were not available when this manuscript file was generated."
    else:
        row = dhx29_meta.iloc[0]
        cross_dataset_text = (
            "A pilot cross-dataset expansion added two Myxoma/poxvirus effector contrasts and found DHX29 "
            f"in the same positive direction across all three pilot contexts, with significance in "
            f"{int(row['significant_dataset_count'])} of {int(row['dataset_count'])} datasets. "
            f"The mean DHX29 log2 fold change across pilot contexts was {row['mean_log2_fold_change']:.3f}, "
            f"with direction consistency {row['direction_consistency']:.3f}."
        )
    if mech_summary.empty or mech_priority.empty:
        mechanistic_text = "Mechanistic prioritization outputs were not available when this manuscript file was generated."
        top_candidate_text = "No mechanistic priority table was available."
    else:
        metrics = dict(zip(mech_summary["metric"], mech_summary["value"], strict=False))
        dhx29_priority = mech_priority[mech_priority["gene_symbol"].eq("DHX29")].iloc[0]
        top_candidate = mech_priority.iloc[0]
        top_enrichment = mech_enrichment.iloc[0] if not mech_enrichment.empty else None
        enrichment_sentence = ""
        if top_enrichment is not None:
            enrichment_sentence = (
                f" The strongest curated module enrichment was {top_enrichment['query_set']} / "
                f"{top_enrichment['module']} (overlap {int(top_enrichment['overlap_size'])}; "
                f"adjusted p-value {top_enrichment['padj']:.3g})."
            )
        mechanistic_text = (
            f"GENCODE v44 biotype annotation reduced the result set to "
            f"{int(float(metrics.get('curated_host_protein_coding_translation_targets', 0)))} curated host protein-coding "
            f"translation targets, of which {int(float(metrics.get('curated_significant_with_abs_log2fc_ge_threshold', 0)))} "
            f"were significant at adjusted p-value < 0.05 with absolute log2 fold change at least 0.5."
            f"{enrichment_sentence} DHX29 ranked {int(dhx29_priority['priority_rank'])} overall and "
            f"{int(float(metrics.get('dhx29_helicase_priority_rank', 0)))} within the DHX/DDX helicase subset; "
            f"within infected samples, DHX29 showed a positive exploratory association with viral RNA burden "
            f"(Pearson r={dhx29_priority['pearson_r_viral_load']:.3f}; p={dhx29_priority['pearson_p_viral_load']:.3g}; n={int(dhx29_priority['n_infected_samples'])})."
        )
        top_candidate_text = (
            f"The highest composite mechanistic-priority candidate was {top_candidate['gene_symbol']} "
            f"(score {top_candidate['candidate_priority_score']:.3f}; log2 fold change "
            f"{top_candidate['log2FoldChange']:.3f}; adjusted p-value {top_candidate['padj']:.3g})."
        )

    manuscript = rf"""\documentclass[11pt]{{article}}
\usepackage[margin=1in]{{geometry}}
\usepackage{{graphicx}}
\usepackage{{booktabs}}
\usepackage{{hyperref}}
\usepackage{{longtable}}
\usepackage{{natbib}}
\usepackage{{setspace}}
\doublespacing

\title{{PoxHostAtlas: A reproducible cross-study atlas of host gene-expression remodeling during poxvirus infection}}
\author{{Morsal Mubarak\\Justice High School, Falls Church, VA \and Deon Nguyen\\Arizona State University, Tempe, AZ}}
\date{{\today}}

\begin{{document}}
\maketitle

\begin{{abstract}}
Poxviruses maintain viral protein synthesis while suppressing host translation, but the broader transcriptomic behavior of host translation-initiation factors remains an interpretable secondary-analysis target in recently published sequencing datasets. Here, I analyze public RNA-sequencing count data from GSE278320, generated by Park and colleagues, using a reproducible Python and PyDESeq2 workflow \citep{{park2025poxvirus,love2014deseq2}}. The analysis focuses on parental HAP1 total-fraction samples comparing Vaccinia virus infection with mock controls, then filters differential-expression results for DHX, EIF, RPS, and RPL genes. After low-count filtering, the model evaluated 21,835 genes and identified 184 statistically significant translation-factor impacts at an adjusted p-value threshold of 0.05. DHX29 was significantly upregulated in the primary model (log2 fold change {dhx29['log2FoldChange']:.3f}; adjusted p-value {dhx29['padj']:.3g}), while multiple eIF3-family genes showed significant negative fold changes. A strict-title sensitivity model excluding two GEO samples with title/filename fraction disagreement retained the DHX29 direction of effect (log2 fold change {strict_dhx29['log2FoldChange_strict']:.3f}). A pilot expansion across two additional public Myxoma/poxvirus effector contrasts tested whether the DHX29 signal remained directionally consistent outside the anchor Vaccinia dataset. A mechanistic evidence layer then used GENCODE v44 host biotypes, viral RNA burden, curated module enrichment, and STRING physical-network evidence to prioritize protein-coding host factors \citep{{frankish2021gencode,szklarczyk2023string}}. This work is framed as a secondary computational analysis, not as the primary discovery of the GSE278320 dataset or the RACK1/eIF3 mechanism reported by Park et al. Its contribution is an auditable pipeline and a DHX-family-focused prioritization layer for mentor-guided biological follow-up.
\end{{abstract}}

\section{{Introduction}}
Vaccinia virus infection causes substantial remodeling of host translation. Park and colleagues recently reported that distinct non-canonical translation-initiation modes arise for selected host and viral mRNAs during poxvirus-induced shutoff, emphasizing RACK1, eIF3, JUN, polysome occupancy, and cryo-electron microscopy evidence for virus-modified ribosomal states \citep{{park2025poxvirus}}. Because that study generated a rich public count dataset, it also creates an opportunity for transparent secondary analysis focused on questions outside the original center of emphasis.

The present study asks a narrower computational question: what does the public GSE278320 count dataset reveal about host DHX-family and broader translation-factor transcriptomic disruption during Vaccinia virus infection? This focus is motivated by the biological relevance of RNA helicases such as DHX29 to translation initiation and by the value of generating mentor-reviewable hypotheses for poxvirus-host translation biology. The analysis deliberately avoids claiming primary discovery of the dataset or the RACK1/eIF3 axis.

\section{{Methods}}
\subsection{{Dataset and Sample Selection}}
The analysis used public NCBI GEO accession GSE278320. The primary contrast was parental HAP1 total-fraction RNA-seq, Vaccinia virus versus mock. Official supplementary per-sample gene-count files were downloaded directly from NCBI GEO. Supplementary count filenames were used as the source of truth for fraction selection because two samples, GSM8544828 and GSM8544830, contained GEO title/filename fraction disagreement. These conflicts were retained in the sample manifest rather than silently corrected.

\subsection{{Differential Expression}}
Raw gene-count files were merged into a gene-row count matrix, filtered to remove genes with fewer than 10 total reads across selected samples, transposed to the samples-by-genes orientation required by PyDESeq2, and modeled with design \texttt{{\string~infection}}. The reported contrast is Vaccinia virus versus mock. Adjusted p-values were taken from the PyDESeq2 results table.

\subsection{{Targeted Translation-Factor Analysis}}
The full differential-expression table was filtered for gene symbols beginning with DHX, EIF, RPS, or RPL. Significant translation-factor impacts were defined as adjusted p-value less than 0.05. Publication figures were generated as PNG, PDF, and SVG files. A strict-title sensitivity model excluded GSM8544828 and GSM8544830 to test whether the metadata conflict changed core interpretation.

\subsection{{Mechanistic Prioritization}}
To address the biological interpretability limits of prefix-based target filtering, host genes were annotated with GENCODE v44 biotypes and a curated host-only protein-coding target set was generated \citep{{frankish2021gencode}}. Viral genes were retained in full audit tables but excluded from the host-only curated volcano plot. The mechanistic layer added four evidence classes: primary differential expression, strict-title sensitivity, cross-dataset directionality, and exploratory correlation between host transcript abundance and summed normalized Vaccinia viral RNA burden across infected samples. Curated translation modules were tested by hypergeometric enrichment, and a high-confidence STRING physical-interaction subnetwork was retrieved for top prioritized host protein-coding candidates \citep{{szklarczyk2023string}}.

\section{{Results}}
\subsection{{Primary Model Performance}}
The primary model included 9 samples, with 3 mock and 6 Vaccinia virus samples. After low-count filtering, 21,835 genes were modeled. The targeted translation-factor screen identified 465 DHX, EIF, RPS, or RPL genes, of which 184 met the adjusted p-value threshold of 0.05.

\subsection{{DHX29 and eIF-Family Signals}}
DHX29 was significantly upregulated in the primary model, with log2 fold change {dhx29['log2FoldChange']:.3f} and adjusted p-value {dhx29['padj']:.3g}. Several eIF3-family genes showed significant negative fold changes, including:
\begin{{center}}
\begin{{tabular}}{{lrr}}
\toprule
Gene & log2 fold change & adjusted p-value \\
\midrule
"""
    for _, row in top_eif3.iterrows():
        manuscript += f"{row['gene_symbol']} & {row['log2FoldChange']:.3f} & {row['padj']:.3g} \\\\\n"
    manuscript += rf"""\bottomrule
\end{{tabular}}
\end{{center}}

\subsection{{Metadata Sensitivity}}
The strict-title sensitivity model excluded GSM8544828 and GSM8544830, leaving 2 mock and 5 Vaccinia virus samples. DHX29 retained a positive direction of effect in this sensitivity model, with strict-title log2 fold change {strict_dhx29['log2FoldChange_strict']:.3f} and adjusted p-value {strict_dhx29['padj_strict']:.3g}. This supports using the primary model as a hypothesis-generating analysis while transparently documenting the metadata issue.

\subsection{{Pilot Cross-Dataset Expansion}}
{cross_dataset_text} The pilot expansion is included to test reproducibility of the DHX/helicase-centered signal, not to make a final universality claim. The complete pilot registry and meta-signature table are included as supplementary release tables.

\subsection{{Host-Curated Mechanistic Evidence}}
{mechanistic_text} {top_candidate_text} These results strengthen the analysis by separating host protein-coding targets from viral genes and pseudogene-like entries, then ranking candidates through convergent statistical, sensitivity, cross-dataset, viral-load, and protein-network evidence layers.

\section{{Discussion}}
This secondary analysis recovers a coherent translation-factor disruption pattern from public GSE278320 count data and prioritizes DHX-family signals for biological follow-up. The DHX29 result is useful as a mentor-review point because it connects a public Vaccinia virus dataset to a DHX-helicase-centered research question and to prior evidence that host RNA helicases, including DHX29, can regulate Myxoma virus replication phenotypes \citep{{myxoma2017dhx}}. The proposed biological model is that poxvirus infection produces a coordinated host translation-remodeling state in which eIF3/eIF4/ribosomal changes coexist with a DHX/DDX helicase response. DHX29 is a plausible focal candidate because DHX29 and eIF3 cooperate during scanning on structured mRNAs, giving this transcriptomic signal a mechanistic entry point for follow-up \citep{{pisareva2016dhx29}}. At the same time, the analysis remains transcript-level and cannot by itself establish protein abundance, ribosome occupancy, or direct mechanism. The strongest near-term interpretation is therefore not that DHX29 is mechanistically proven to mediate poxvirus translation, but that DHX29 and associated translation factors are high-priority candidates for comparison with existing poxvirus and Myxoma virus hypotheses.

\section{{Limitations}}
The primary model uses total-fraction count files and does not directly model ribosome profiling, polysome occupancy, or protein abundance. Two GEO records contain title/filename fraction disagreement, handled by filename-based selection and strict-title sensitivity analysis. The strict-title model has fewer mock samples and should be treated as sensitivity evidence rather than a replacement for the primary model. Viral-load correlations are exploratory because only six infected samples are available in the primary contrast. STRING interactions support network context but do not prove condition-specific physical binding. Additional wet-lab or orthogonal computational validation would be required for mechanistic claims.

\section{{Data and Code Availability}}
All code, selected outputs, validation reports, hash receipts, manuscript files, and release assets are included in the accompanying reproducible research package. Raw count files are obtained from NCBI GEO accessions GSE278320, GSE287860, and GSE288000. The original primary publication for the anchor dataset is Park et al. \citep{{park2025poxvirus}}.

\bibliographystyle{{plainnat}}
\bibliography{{references}}
\end{{document}}
"""

    references = """@article{park2025poxvirus,
  title = {Distinct non-canonical translation initiation modes arise for specific host and viral mRNAs during poxvirus-induced shutoff},
  author = {Park, Chorong and Ferrell, Aaron J. and Meade, Nathan and Shen, Peter S. and Walsh, Derek and others},
  journal = {Nature Microbiology},
  volume = {10},
  pages = {1535--1549},
  year = {2025},
  doi = {10.1038/s41564-025-02009-4}
}

@article{love2014deseq2,
  title = {Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2},
  author = {Love, Michael I. and Huber, Wolfgang and Anders, Simon},
  journal = {Genome Biology},
  volume = {15},
  pages = {550},
  year = {2014},
  doi = {10.1186/s13059-014-0550-8}
}

@article{frankish2021gencode,
  title = {GENCODE 2021},
  author = {Frankish, Adam and Diekhans, Mark and Ferreira, Anne-Maud and Johnson, Rory and Jungreis, Irwin and Loveland, Jane and Mudge, Jonathan M. and Sisu, Cristina and Wright, James and Armstrong, Joel and others},
  journal = {Nucleic Acids Research},
  volume = {49},
  number = {D1},
  pages = {D916--D923},
  year = {2021},
  doi = {10.1093/nar/gkaa1087}
}

@article{szklarczyk2023string,
  title = {The STRING database in 2023: protein-protein association networks and functional enrichment analyses for any sequenced genome of interest},
  author = {Szklarczyk, Damian and Kirsch, Rebecca and Koutrouli, Mikaela and Nastou, Katerina and Mehryary, Farrokh and Hachilif, Rebecca and Gable, Annika L. and Fang, Tao and Doncheva, Nadezhda T. and Pyysalo, Sampo and others},
  journal = {Nucleic Acids Research},
  volume = {51},
  number = {D1},
  pages = {D638--D646},
  year = {2023},
  doi = {10.1093/nar/gkac1000}
}

@article{pisareva2016dhx29,
  title = {DHX29 and eIF3 cooperate in ribosomal scanning on structured mRNAs during translation initiation},
  author = {Pisareva, Vera P. and Pisarev, Andrey V.},
  journal = {RNA},
  volume = {22},
  number = {12},
  pages = {1859--1870},
  year = {2016},
  doi = {10.1261/rna.057851.116}
}

@article{myxoma2017dhx,
  title = {Identification of host DEAD-box RNA helicases that regulate cellular tropism of oncolytic Myxoma virus in human cancer cells},
  author = {Rahman, Masmudur M. and others},
  journal = {Scientific Reports},
  volume = {7},
  pages = {15710},
  year = {2017},
  doi = {10.1038/s41598-017-15941-1}
}
"""

    (RELEASE_DIR / "manuscript" / "main.tex").write_text(manuscript, encoding="utf-8")
    (RELEASE_DIR / "manuscript" / "references.bib").write_text(references, encoding="utf-8")


def copy_release_assets() -> None:
    for script in sorted((REPO_ROOT / "scripts").glob("*.py")):
        copy_file(script, RELEASE_DIR / "code" / "scripts" / script.name)
    for filename in ["requirements.txt", "requirements-lock.txt", "README.md"]:
        copy_file(REPO_ROOT / filename, RELEASE_DIR / "code" / filename)

    for source in sorted((RAW_DIR / "count_files").glob("*.txt.gz")):
        copy_file(source, RELEASE_DIR / "data" / "raw" / "count_files" / source.name)
    for source in sorted((RAW_DIR / "geo_metadata").glob("*")):
        if source.is_file():
            copy_file(source, RELEASE_DIR / "data" / "raw" / "geo_metadata" / source.name)
    for source in sorted(EXTERNAL_RAW_DIR.glob("*.gz")):
        copy_file(source, RELEASE_DIR / "data" / "external" / "expansion_raw" / source.name)

    for source in sorted(PROCESSED_DIR.glob("*.csv")):
        copy_file(source, RELEASE_DIR / "data" / "processed" / source.name)

    figure_map = {
        "volcano_plot": "Figure_2_volcano_plot",
        "translation_heatmap": "Figure_3_translation_heatmap",
        "poxvirus_translation_factor_meta_signature": "Figure_4_cross_dataset_meta_signature",
        "host_only_curated_volcano": "Figure_5_host_only_curated_volcano",
        "candidate_evidence_matrix": "Figure_6_mechanistic_evidence_matrix",
        "dhx_ddx_helicase_priority": "Figure_7_dhx_ddx_helicase_priority",
        "dhx29_viral_load_correlation": "Figure_8_dhx29_viral_load_correlation",
        "module_enrichment_barplot": "Figure_9_module_enrichment",
        "string_priority_network": "Figure_10_string_priority_network",
    }
    meta_figure_map = {
        "meta_volcano_pan_poxvirus": "Figure_11_meta_volcano_pan_poxvirus",
        "forest_DHX29": "Figure_12_forest_DHX29",
        "evidence_heatmap_top_tf": "Figure_13_evidence_heatmap_top_tf",
        "forest_DHX15": "Figure_14_forest_DHX15",
        "forest_EIF4B": "Figure_15_forest_EIF4B",
    }
    synthesis_figure_map = {
        "host_background_robustness_heatmap": "Figure_16_host_background_robustness",
        "vaccinia_temporal_dynamics": "Figure_17_vaccinia_temporal_dynamics",
        "module_concordance_summary": "Figure_18_module_concordance_summary",
        "pipeline_overview": "Figure_19_pipeline_overview",
        "cross_study_rank_correlation": "Figure_20_cross_study_rank_correlation",
        "extended_evidence_heatmap": "Figure_21_extended_evidence_heatmap",
    }
    for stem, release_stem in figure_map.items():
        for suffix in ("png", "pdf", "svg"):
            if stem.startswith("poxvirus_"):
                source_dir = META_DIR
            elif stem in {
                "host_only_curated_volcano",
                "candidate_evidence_matrix",
                "dhx_ddx_helicase_priority",
                "dhx29_viral_load_correlation",
                "module_enrichment_barplot",
                "string_priority_network",
            }:
                source_dir = MECH_DIR
            else:
                source_dir = FIGURE_DIR
            source = source_dir / f"{stem}.{suffix}"
            if source.exists():
                copy_file(source, RELEASE_DIR / "figures" / f"{release_stem}.{suffix}")

    ATLAS_DIR = RESULTS_DIR / "figures" / "atlas"
    for i, stem in [(1, "Figure1_atlas"), (2, "Figure2_meta"), (3, "Figure3_modules"),
                    (4, "Figure4_ml"), (5, "Figure5_network"), (6, "Figure6_validation")]:
        for suffix in ("png", "pdf"):
            source = ATLAS_DIR / f"{stem}.{suffix}"
            if source.exists():
                copy_file(source, RELEASE_DIR / "figures" / f"Main_Figure_{i}_{stem.split('_',1)[1]}.{suffix}")

    META_ANALYSIS_DIR = RESULTS_DIR / "meta_analysis" / "figures"
    SYNTHESIS_DIR = RESULTS_DIR / "synthesis" / "figures"
    for stem, release_stem in meta_figure_map.items():
        for suffix in ("png", "pdf", "svg"):
            source = META_ANALYSIS_DIR / f"{stem}.{suffix}"
            if source.exists():
                copy_file(source, RELEASE_DIR / "figures" / f"{release_stem}.{suffix}")
    for stem, release_stem in synthesis_figure_map.items():
        for suffix in ("png", "pdf", "svg"):
            source = SYNTHESIS_DIR / f"{stem}.{suffix}"
            if source.exists():
                copy_file(source, RELEASE_DIR / "figures" / f"{release_stem}.{suffix}")

    manuscript_src = REPO_ROOT / "docs" / "manuscript" / "main.tex"
    if manuscript_src.exists():
        copy_file(manuscript_src, RELEASE_DIR / "manuscript" / "main.tex")

    table_sources = {
        PROCESSED_DIR / "sample_manifest.csv": "Table_1_sample_manifest.csv",
        TABLE_DIR / "top_downregulated_translation_factors.csv": "Table_2_top_downregulated_translation_factors.csv",
        TABLE_DIR / "top_upregulated_translation_factors.csv": "Table_3_top_upregulated_translation_factors.csv",
        TABLE_DIR / "dge_results_full.csv": "Supplementary_Table_1_full_dge_results.csv",
        TABLE_DIR / "translation_factors_all.csv": "Supplementary_Table_2_translation_factors_all.csv",
        TABLE_DIR / "translation_factors_impact.csv": "Supplementary_Table_3_translation_factors_impact.csv",
        SENSITIVITY_DIR / "sensitivity_primary_vs_strict_translation.csv": "Supplementary_Table_4_primary_vs_strict_sensitivity.csv",
        SENSITIVITY_DIR / "translation_factors_impact_strict_titles.csv": "Supplementary_Table_5_strict_title_translation_factors_impact.csv",
        META_DIR / "poxvirus_translation_factor_meta_signature.csv": "Supplementary_Table_6_cross_dataset_meta_signature.csv",
        META_DIR / "GSE287860_translation_targets.csv": "Supplementary_Table_7_GSE287860_translation_targets.csv",
        META_DIR / "GSE288000_NTC_translation_targets.csv": "Supplementary_Table_8_GSE288000_NTC_translation_targets.csv",
        META_DIR / "cross_dataset_pilot_registry.csv": "Supplementary_Table_9_cross_dataset_pilot_registry.csv",
        MECH_DIR / "mechanistic_priority_rank.csv": "Supplementary_Table_10_mechanistic_priority_rank.csv",
        MECH_DIR / "top_mechanistic_candidates.csv": "Supplementary_Table_11_top_mechanistic_candidates.csv",
        MECH_DIR / "dhx_ddx_helicase_priority_rank.csv": "Supplementary_Table_12_dhx_ddx_helicase_priority_rank.csv",
        MECH_DIR / "curated_module_enrichment.csv": "Supplementary_Table_13_curated_module_enrichment.csv",
        MECH_DIR / "host_gene_viral_load_correlations.csv": "Supplementary_Table_14_host_gene_viral_load_correlations.csv",
        MECH_DIR / "viral_load_by_sample.csv": "Supplementary_Table_15_viral_load_by_sample.csv",
        MECH_DIR / "string_network_edges.csv": "Supplementary_Table_16_string_network_edges.csv",
        MECH_DIR / "string_network_nodes.csv": "Supplementary_Table_17_string_network_nodes.csv",
        MECH_DIR / "translation_targets_biotype_audit.csv": "Supplementary_Table_18_translation_targets_biotype_audit.csv",
        MECH_DIR / "mechanistic_summary_metrics.csv": "Supplementary_Table_19_mechanistic_summary_metrics.csv",
        RESULTS_DIR / "meta_analysis" / "meta_summary.csv": "Supplementary_Table_20_meta_summary.csv",
        RESULTS_DIR / "meta_analysis" / "top_conserved_translation_factors.csv": "Supplementary_Table_21_top_conserved_translation_factors.csv",
        RESULTS_DIR / "expanded" / "expanded_contrast_registry.csv": "Supplementary_Table_22_expanded_contrast_registry.csv",
        RESULTS_DIR / "expanded" / "GSE288000_host_background_robustness.csv": "Supplementary_Table_23_host_background_robustness.csv",
        RESULTS_DIR / "synthesis" / "module_concordance_summary.csv": "Supplementary_Table_24_module_concordance_summary.csv",
        RESULTS_DIR / "synthesis" / "composite_evidence_score.csv": "Supplementary_Table_25_composite_evidence_score.csv",
        RESULTS_DIR / "synthesis" / "novel_findings_summary.csv": "Supplementary_Table_26_novel_findings_summary.csv",
        RESULTS_DIR / "synthesis" / "cross_study_spearman_correlation.csv": "Supplementary_Table_27_cross_study_spearman_correlation.csv",
        RESULTS_DIR / "meta_analysis" / "leave_one_study_out.csv": "Supplementary_Table_28_leave_one_study_out.csv",
        RESULTS_DIR / "meta_analysis" / "heterogeneity_classification.csv": "Supplementary_Table_29_heterogeneity_classification.csv",
        RESULTS_DIR / "meta_analysis" / "rank_concordance_matrix.csv": "Supplementary_Table_30_rank_concordance_matrix.csv",
        RESULTS_DIR / "ml" / "leave_dataset_out_performance.csv": "Supplementary_Table_31_leave_dataset_out_performance.csv",
        RESULTS_DIR / "ml" / "leave_virus_out_performance.csv": "Supplementary_Table_32_leave_virus_out_performance.csv",
        RESULTS_DIR / "ml" / "feature_importance_consensus.csv": "Supplementary_Table_33_feature_importance_consensus.csv",
        RESULTS_DIR / "ml" / "ablation_results.csv": "Supplementary_Table_34_ablation_results.csv",
        RESULTS_DIR / "ml" / "negative_control_results.csv": "Supplementary_Table_35_negative_control_results.csv",
        RESULTS_DIR / "network" / "network_node_scores.csv": "Supplementary_Table_36_network_node_scores.csv",
        RESULTS_DIR / "synthesis" / "final_candidate_ranking.csv": "Supplementary_Table_37_final_candidate_ranking.csv",
        DOCS_DATASET_REGISTRY: "Supplementary_Table_38_dataset_registry.csv",
    }
    for source, name in table_sources.items():
        if source.exists():
            copy_file(source, RELEASE_DIR / "tables" / name)

    copy_file(REPORT_DIR / "validation_report.md", RELEASE_DIR / "provenance" / "validation_report.md")
    copy_file(REPORT_DIR / "file_hashes_sha256.csv", RELEASE_DIR / "provenance" / "analysis_file_hashes_sha256.csv")
    copy_file(REPO_ROOT / "docs" / "SOURCES.md", RELEASE_DIR / "provenance" / "source_urls.md")
    copy_file(REPO_ROOT / "docs" / "NOVELTY_AUDIT.md", RELEASE_DIR / "provenance" / "novelty_audit.md")
    copy_file(REPO_ROOT / "docs" / "EXPANSION_DATASET_REGISTRY.csv", RELEASE_DIR / "provenance" / "expansion_dataset_registry.csv")
    copy_file(REPO_ROOT / "docs" / "public_dataset_discovery_results.csv", RELEASE_DIR / "provenance" / "public_dataset_discovery_results.csv")
    copy_file(REPO_ROOT / "requirements-lock.txt", RELEASE_DIR / "provenance" / "environment_lock.txt")


def write_release_readme() -> None:
    text = """# PoxHostAtlas - Reproducible Release Package

**PoxHostAtlas** is a pan-poxvirus computational atlas and interpretable predictive model of conserved host translation/RNA-helicase remodeling. This package contains the manuscript, six main figures, all DGE/meta-analysis/ML/network tables, the systematic dataset registry, validation report, and SHA-256 provenance.

Core claim: across a systematically curated atlas of public poxvirus host-response datasets, infection reproducibly remodels host translation/helicase programs; this signature predicts infection in unseen studies, survives ablation and falsifiable negative controls, and prioritizes **DHX15** (conserved, Tier 1) and **DHX29** (context-dependent, Class II) for experimental testing.

## Reproduce

```bash
bash run_all.sh
```

## Methods Stack

- Systematic dataset discovery + A/B/C tiering (`docs/dataset_registry.csv`, `docs/dataset_inclusion_criteria.md`).
- Harmonized PyDESeq2 differential expression across 10 contrasts (4 Tier-A datasets: 3 raw-count reanalyzed + 1 author-DE-table time course).
- DerSimonian-Laird random-effects meta-analysis + leave-one-study-out robustness + Class I-IV heterogeneity classification + signed-rank concordance (Kendall's W).
- Leave-one-DATASET-out interpretable ML (elastic-net, RF, gradient boosting, linear SVM) with bootstrap 95% CIs, expression/variance-matched null gene sets, ablation with permutation significance, and negative controls.
- STRING network centrality + multi-evidence Final Evidence Score; external GSE185520 validation.

## Key Results

- Pan-poxvirus meta-analysis: 65 translation factors at meta-FDR < 0.05 (DHX/DDX up; eIF3/eIF4B/ribosomal down).
- DHX15: top conserved helicase (meta-FDR 3.6e-19; LOSO-robust; Tier 1). DHX29: context-dependent (meta-FDR ~0.055; Class II, I^2 ~ 80%; Tier 3).
- LODO prediction (linear SVM, translation factors): balanced accuracy 0.79 (95% CI 0.67-0.90); ROC-AUC 0.88 (0.76-0.98).
- Ablation: removing top-25 conserved genes degrades prediction more than 99.8% of 1,000 matched random ablations (p=0.003).
- Negative controls: label permutation 0.54 (chance); study-identity 0.93 (batch probe).

## Interpretation Boundary

Computational, reproducible, hypothesis-generating. Sample-level ML uses three raw-count studies (45 samples); leave-virus-out is exploratory (only two virus groups). Does not claim primary dataset discovery (e.g., Park et al. GSE278320), broad poxvirus universality beyond harmonizable Orthopoxvirus/Leporipoxvirus data, or causal proof of DHX15/DHX29 function.
"""
    (RELEASE_DIR / "README_RELEASE.md").write_text(text, encoding="utf-8")


def write_run_log() -> None:
    commands = [
        "python scripts/01_fetch_data.py",
        "python scripts/02_preprocess_counts.py --min-total-count 10",
        "python scripts/03_run_pydeseq2.py --n-cpus 8",
        "python scripts/04_analyze_translation_factors.py --padj 0.05",
        "python scripts/05_generate_figures.py --padj 0.05 --top-n 30",
        "python scripts/06_run_sensitivity_analysis.py --min-total-count 10 --padj 0.05 --n-cpus 8",
        "python scripts/10_run_cross_dataset_pilot.py --min-total-count 10 --padj 0.05 --n-cpus 8",
        "python scripts/12_run_expanded_dge.py --min-total-count 10 --padj 0.05 --n-cpus 8",
        "python scripts/13_run_meta_analysis.py",
        "python scripts/14_run_advanced_synthesis.py",
        "python scripts/16_run_integrative_synthesis.py",
        "python scripts/18_run_robust_meta.py",
        "python scripts/19_run_ml_atlas.py",
        "python scripts/20_run_network_ranking.py",
        "MPLBACKEND=Agg python scripts/21_make_atlas_figures.py",
        "python scripts/07_generate_validation_report.py",
        "python scripts/08_prepare_release.py",
    ]
    text = "# Reproducibility Run Log\n\n" + "\n".join(f"- `{command}`" for command in commands) + "\n"
    (RELEASE_DIR / "provenance" / "run_log.txt").write_text(text, encoding="utf-8")


def write_release_hashes() -> None:
    rows = []
    for path in sorted(RELEASE_DIR.rglob("*")):
        if path.is_file() and path.name != "release_file_hashes_sha256.csv":
            rows.append(
                {
                    "relative_path": path.relative_to(RELEASE_DIR).as_posix(),
                    "size_bytes": path.stat().st_size,
                    "sha256": sha256_file(path),
                }
            )
    pd.DataFrame(rows).to_csv(RELEASE_DIR / "provenance" / "release_file_hashes_sha256.csv", index=False)


def zip_release() -> Path:
    archive_path = REPO_ROOT.parent / f"openbio_sri_full_release_{date.today().isoformat()}.zip"
    if archive_path.exists():
        archive_path.unlink()
    with zipfile.ZipFile(archive_path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(RELEASE_DIR.rglob("*")):
            if path.is_file():
                archive.write(path, Path("release") / path.relative_to(RELEASE_DIR))
    hash_path = archive_path.with_suffix(archive_path.suffix + ".sha256")
    hash_path.write_text(f"{sha256_file(archive_path)}  {archive_path.name}\n", encoding="utf-8")
    return archive_path


def main() -> None:
    reset_release_tree()
    copy_release_assets()
    write_pipeline_schematic()
    (RELEASE_DIR / "manuscript").mkdir(parents=True, exist_ok=True)
    for _pdf in ("PoxHostAtlas_manuscript.pdf", "PoxHostAtlas_supplementary.pdf"):
        _src = REPO_ROOT / "docs" / "manuscript" / _pdf
        if _src.exists():
            copy_file(_src, RELEASE_DIR / "manuscript" / _pdf)
    write_release_readme()
    write_run_log()
    write_release_hashes()
    archive_path = zip_release()
    print(f"Wrote release tree: {RELEASE_DIR}")
    print(f"Wrote release archive: {archive_path}")
    print(f"Wrote archive hash: {archive_path.with_suffix(archive_path.suffix + '.sha256')}")


if __name__ == "__main__":
    main()
