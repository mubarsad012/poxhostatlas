# PoxHostAtlas

**A reproducible cross-study atlas of host gene-expression remodeling during poxvirus infection**

- Author: Morsal Mubarak (Justice High School, Falls Church, VA) · [morsalmubarak9@gmail.com](mailto:morsalmubarak9@gmail.com)
- Mentor / last author: Deon Nguyen (Arizona State University, Tempe, AZ)
- OpenBio Student Research Institute · computational secondary analysis
- Full write-up: `docs/manuscript/PoxHostAtlas_manuscript.pdf` · `docs/manuscript/PoxHostAtlas_supplementary.pdf`

## Question

- Genome-wide, unbiased, family-agnostic
- Which host genes/pathways are reproducibly remodeled across independent poxvirus datasets?
- Which host hubs are most consistently targeted?
- RNA-helicase remodeling = one pre-specified sub-question, not the headline

## Key findings

- **1,834** host genes reproducibly remodeled (871 up, 963 down)
- Reproducible = ≥3 studies · concordant direction · meta-FDR < 0.05 · |log2FC| ≥ 0.5
- Two dominant programs:
  - ↑ inflammatory / NF-kB signaling (5.4×, FDR 3e-7) — IL6, IL1B, CXCL8, BIRC3
  - ↓ cell-cycle + ribosomal-protein genes (5.7×, FDR 3e-7) = host shutoff
- Vaccinia-only sensitivity: shutoff = cleanest signal; inflammation real but effector-amplified
- Helicases = small module (29/1,834 = 1.6%); **DHX15** top helicase (meta-FDR 3.6e-19, top 0.6%); DHX29 context-dependent
- Central hubs: IL6, PTPRC, MAFB
- Portable to held-out studies (balanced accuracy 0.79; permuted 0.54)
- External validation (Matía 2024): single-cell shutoff in 2 cell types + transcript↔protein concordance (ρ = 0.31)

## Datasets


| Tier | Datasets                                   | Role                                    |
| ---- | ------------------------------------------ | --------------------------------------- |
| A    | GSE278320, GSE284044, GSE287860, GSE288000 | quantitative integration (10 contrasts) |
| B    | GSE185520                                  | external validation                     |
| C    | 7 (long-read viral, expansion targets)     | context only                            |


- 12 screened · 16,841 host genes in ≥3 studies
- Registry: `docs/dataset_registry.csv` · criteria: `docs/dataset_inclusion_criteria.md` · excluded: `docs/excluded_datasets.csv`

## Pipeline

- Discovery + pre-registered tiering
- Harmonized differential expression — PyDESeq2; Vero→human orthologues (GENCODE v44)
- Genome-wide DerSimonian–Laird random-effects meta-analysis — I², concordance, meta-FDR, leave-one-study-out
- Host-program over-representation (hypergeometric) + vaccinia-only layer
- Family-agnostic cross-study co-expression network + centrality + communities
- Leave-one-dataset-out portability stress-test — permutation, ablation, negative controls
- Independent external validation — single-cell + proteomics

## Modules (script → output)

- `22_genomewide_reproducibility.py` → `results/genomewide/` — reproducible host genes, program enrichment, vaccinia-only
- `23_genomewide_network_ranking.py` → `results/genomewide/network/` — hub ranking, communities
- `18_run_robust_meta.py` — leave-one-study-out, Class I–IV heterogeneity, Spearman + Kendall's W
- `13_run_meta_analysis.py`, `20_run_network_ranking.py` — DHX15/DHX29 forests, STRING, multi-evidence rank (sub-analysis)
- `19_run_ml_atlas.py` — portability stress-test (secondary; not a prediction product)
- `21_make_atlas_figures.py` → `results/genomewide/figures/` — atlas, volcano, enrichment, network, forests
- `12_run_expanded_dge.py`, `14_run_advanced_synthesis.py`, `16_run_integrative_synthesis.py` → `results/expanded/`, `results/synthesis/`

## Setup

```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

## Run

```bash
bash run_all.sh          # full release run
```

```bash
python scripts/01_fetch_data.py
python scripts/02_preprocess_counts.py --min-total-count 10
python scripts/03_run_pydeseq2.py --n-cpus 8
python scripts/22_genomewide_reproducibility.py
python scripts/23_genomewide_network_ranking.py
python scripts/07_generate_validation_report.py
python scripts/08_prepare_release.py
```

## Key outputs

- `results/genomewide/` — ranked reproducible genes, `top50_reproducible_{up,down}.csv`, `program_enrichment.csv`, `program_direction_summary{,_vaccinia_only}.csv`, `helicase_rank_in_context.csv`
- `results/genomewide/network/` — `genomewide_network_hub_ranking.csv`, `top_host_hubs.csv`, `network_communities.csv`
- `results/genomewide/figures/` — `genomewide_volcano.png`, `program_enrichment.png`, network + forest plots
- `results/reports/validation_report.md` — QC + SHA-256 manifest
- `docs/manuscript/` — final manuscript + supplementary (PDF)
- `docs/OPENBIO_SRI_SUBMISSION_BRIEF.md` — mentor/SRI brief
- `release/` — packaged release (data, figures, tables, provenance, hashes, code)

## Scope

- Computational · reproducible · hypothesis-generating
- No claim of primary dataset discovery, the RACK1/eIF3 mechanism, or causal host-gene function
- Transcript-level signals + one independent proteomic validation

## Data availability

- Public accessions: GSE278320, GSE284044, GSE287860, GSE288000, GSE185520
- Code, registry, tables, figures released with SHA-256 manifest

## References

1. Park C, et al. Nat Microbiol. 2025;10:1535-1549.
2. Rahman MM, et al. Sci Rep. 2017;7:15710.
3. Pisareva VP, Pisarev AV. RNA. 2016;22(12):1859-1870.
4. Tombácz D, et al. Pathogens. 2021;10(8):919.
5. Tombácz D, et al. GigaScience. 2018;7(12):giy139.
6. Yang Z, et al. J Virol. 2015;89(13):6874-6886.
7. Matía A, et al. bioRxiv. 2024. doi:10.1101/2024.01.13.575413.
8. Muzellec B, et al. Bioinformatics. 2023;39(9):btad547.
9. Love MI, et al. Genome Biol. 2014;15:550.
10. Szklarczyk D, et al. Nucleic Acids Res. 2023;51(D1):D638-D646.
11. Frankish A, et al. Nucleic Acids Res. 2021;49(D1):D916-D923.

