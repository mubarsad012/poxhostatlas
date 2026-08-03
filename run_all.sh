#!/usr/bin/env bash
set -euo pipefail

PYTHON_BIN="${PYTHON:-python}"
if [[ -x "venv/bin/python" ]]; then
  PYTHON_BIN="venv/bin/python"
fi

"${PYTHON_BIN}" -m py_compile scripts/*.py
"${PYTHON_BIN}" scripts/01_fetch_data.py
"${PYTHON_BIN}" scripts/02_preprocess_counts.py --min-total-count 10
"${PYTHON_BIN}" scripts/03_run_pydeseq2.py --n-cpus 8
"${PYTHON_BIN}" scripts/04_analyze_translation_factors.py --padj 0.05
"${PYTHON_BIN}" scripts/05_generate_figures.py --padj 0.05 --top-n 30
"${PYTHON_BIN}" scripts/06_run_sensitivity_analysis.py --min-total-count 10 --padj 0.05 --n-cpus 8
"${PYTHON_BIN}" scripts/10_run_cross_dataset_pilot.py --min-total-count 10 --padj 0.05 --n-cpus 8
"${PYTHON_BIN}" scripts/11_run_mechanistic_prioritization.py --padj 0.05 --effect-size 0.5 --top-n 30
"${PYTHON_BIN}" scripts/12_run_expanded_dge.py --min-total-count 10 --padj 0.05 --n-cpus 8
"${PYTHON_BIN}" scripts/13_run_meta_analysis.py
MPLBACKEND=Agg "${PYTHON_BIN}" scripts/22_genomewide_reproducibility.py
MPLBACKEND=Agg "${PYTHON_BIN}" scripts/23_genomewide_network_ranking.py
MPLBACKEND=Agg "${PYTHON_BIN}" scripts/24_matia_external_validation.py
"${PYTHON_BIN}" scripts/14_run_advanced_synthesis.py
"${PYTHON_BIN}" scripts/16_run_integrative_synthesis.py
"${PYTHON_BIN}" scripts/18_run_robust_meta.py
"${PYTHON_BIN}" scripts/19_run_ml_atlas.py
"${PYTHON_BIN}" scripts/20_run_network_ranking.py
MPLBACKEND=Agg "${PYTHON_BIN}" scripts/21_make_atlas_figures.py
"${PYTHON_BIN}" scripts/07_generate_validation_report.py
"${PYTHON_BIN}" scripts/08_prepare_release.py
