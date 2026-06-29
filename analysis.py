#!/usr/bin/env python3
"""End-to-end reproducible analysis pipeline for PoxHostAtlas.

Executing this script reruns all preceding analysis stages sequentially, replicating every file contained within the results directory.

As the project expands, new stages are added to the pipeline; this version captures the analysis sequence up to the final stage documented here.


Usage:
    python analysis.py                 # run the whole pipeline
    python analysis.py --from 13_run_meta_analysis.py   # resume from a stage
    python analysis.py --list          # list stages and exit

Requirements: see requirements.txt. The data-fetch stage needs internet access.
"""
from __future__ import annotations
import argparse, os, subprocess, sys
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
ROOT = Path(__file__).resolve().parent
SCRIPTS = ROOT / "scripts"

# (label, script, args). Appended to as the project builds.
STAGES = [
    ('Fetch GEO data', '01_fetch_data.py', []),
    ('Preprocess & harmonize counts', '02_preprocess_counts.py', ['--min-total-count', '10']),
    ('Differential expression (PyDESeq2)', '03_run_pydeseq2.py', ['--n-cpus', '4']),
    ('Translation-factor / helicase view', '04_analyze_translation_factors.py', ['--padj', '0.05']),
]


def run_stage(label, script, args):
    print(f"\n=== {label}  ({script}) ===", flush=True)
    subprocess.run([sys.executable, str(SCRIPTS / script), *args], check=True)


def main():
    ap = argparse.ArgumentParser(description="Run the PoxHostAtlas analysis pipeline.")
    ap.add_argument("--from", dest="start", default=None, help="script filename to start from")
    ap.add_argument("--list", action="store_true", help="list stages and exit")
    a = ap.parse_args()
    if a.list:
        for label, script, _ in STAGES:
            print(f"{script:38s} {label}")
        return
    started = a.start is None
    for label, script, args in STAGES:
        started = started or script == a.start
        if started:
            run_stage(label, script, args)


if __name__ == "__main__":
    main()
