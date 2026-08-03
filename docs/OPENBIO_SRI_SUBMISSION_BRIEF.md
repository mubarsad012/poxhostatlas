# OpenBio SRI Submission Brief

- Project: **PoxHostAtlas** — reproducible cross-study atlas of host gene-expression remodeling during poxvirus infection
- Author: Morsal Mubarak (Justice High School, Falls Church, VA) · Mentor / last author: Deon Nguyen (Arizona State University)
- Type: computational secondary analysis · genome-wide · family-agnostic

## Core claim

- Poxvirus infection **reproducibly remodels a broad host program** across independent viruses, hosts, cell types, labs
- ↑ inflammatory / NF-κB signaling · ↓ cell-cycle + ribosomal/translation machinery (host shutoff)
- Helicases = small coherent sub-module; **DHX15** = top reproducible helicase for follow-up
- No pre-selected gene family; no causal-mechanism claim

## Objective

- Turn fragmented single-dataset reanalyses → one reproducible discovery resource
- Question: which host genes/pathways most reproducibly remodeled? which hubs most targeted?
- Pipeline: discovery + pre-registered tiering → harmonized PyDESeq2 DGE → genome-wide DerSimonian–Laird meta-analysis (LOSO + heterogeneity classes) → host-program over-representation → cross-study co-expression hub ranking → vaccinia-only sensitivity → leave-one-dataset-out portability
- Builds on Park et al. (Nat Microbiol 2025) + mentor-flagged transcriptome literature; no primary-discovery claim

## Reframe note (2026-06)

- Earlier version led with RNA helicases + a "predictive model"
- Per mentor feedback: broadened to genome-wide (helicases ≈ 1.6% of signal); ML demoted to secondary stress-test
- Detail: `docs/SUPERVISOR_UPDATE_2026-06-21.md`

## Key results

| Metric | Value |
| --- | --- |
| Datasets screened / Tier A / B / C | 12 / 4 / 1 / 7 |
| Harmonized DGE contrasts | 10 |
| Reproducible host genes (≥3 studies, concordant, meta-FDR<0.05, \|log2FC\|≥0.5) | 1,834 (871 up / 963 down) |
| Translation factors among reproducible genes | 29 / 1,834 (1.6%) |
| Top programs | inflammatory/NF-κB up (5.4×, FDR 3e-7); cell cycle down (5.7×, FDR 3e-7); ribosomal down (24/26) |
| Vaccinia-only (cleanest signal) | ribosome −0.79, cell cycle −0.58, eIF −0.32 down; inflammatory +0.74 |
| DHX15 (top helicase) | pooled log2FC 0.75; meta-FDR 3.6e-19; rank #96/16,841; LOSO-robust |
| DHX29 (context-dependent) | pooled log2FC 0.51; meta-FDR 0.06; high heterogeneity; ~rank 3,400 |
| Co-expression network | 228 nodes / 4,430 edges / 4 communities; hubs IL6, PTPRC, MAFB |
| Portability (leave-one-dataset-out) | balanced acc 0.79; permuted 0.54; study-identity 0.93; ablation p=0.003 |
| External — single-cell, 2 cell types (Matía 2024) | host-read frac 0.998→0.76 (HeLa), 0.999→0.80 (BJ5ta), p≈0 |
| External — proteomics (Matía 2024) | transcript↔protein ρ=0.31 (p=8e-7, n=242); transcript-down→protein-down (−0.21 vs +0.05) |

## Main figures

1. Atlas — discovery (PRISMA-style), scale, virus coverage, pipeline
2. Genome-wide reproducibility — volcano (helicases circled, not pre-selected) + program over-representation
3. Cross-study co-expression network — family-agnostic hub ranking
4. Helicase sub-program — meta-analysis, heterogeneity classes, DHX15/DHX29 forests
5. External validation (Matía 2024) — single-cell shutoff (2 cell types) + transcript↔protein concordance
- Supplementary — portability stress-test (controls, ablation); external GSE185520 validation

## Manuscript

- `docs/manuscript/PoxHostAtlas_manuscript.pdf` — final manuscript
- `docs/manuscript/PoxHostAtlas_supplementary.pdf` — supplementary

## Claim boundary

- Computational · reproducible · hypothesis-generating
- No claim of primary discovery, RACK1/eIF3 mechanism, or causal host-gene function
- Built-in: honest tiering · genome-wide scope · heterogeneity classes · vaccinia-only layer · negative controls
- ML = secondary portability check, not a prediction product

## Edge vs comparators (DROID / CUPNavigator)

- Named, reusable multi-study / virus / host atlas
- Unbiased genome-wide reproducibility map (not one pre-selected family)
- Cross-study generalization via leave-one-dataset-out (not leaky random splits)
- Transparent program enrichment + network hub ranking
- Honest sensitivity analyses · full public release with provenance
