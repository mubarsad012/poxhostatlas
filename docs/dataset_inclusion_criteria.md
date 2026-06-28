# PoxHostAtlas Dataset Criteria For Inclusion / Exclusion

Within this document, the systematic pre-registered criteria being used to build our PoxHostAtlas public poxvirus host-response dataset, will be defined (`docs/dataset_registry.csv`). The goal here is to have an auditable record that is PRISMA-styled that entails how each of the candidate datasets were screened, harmonized, and tiered. 

## Systematic Search

The databases that were queried are as follows: NCBI GEO (`gds`), NCBI SRA, NCBI PubMed (linked supplementary tables), and the GEO FTP supplementary trees. The machine-readable search results are now being archived in `docs/public_dataset_discovery_results.csv`so that we can maintain organization within this project and know where everything is. 

Representative Entrez Direct queries:

```
esearch -db gds -query '(poxvirus OR vaccinia OR myxoma OR monkeypox OR mpox OR \
  orthopoxvirus OR cowpox OR ectromelia OR "orf virus") AND ("RNA-seq" OR \
  transcriptome OR "ribosome profiling" OR polysome)' | efetch -format docsum

esearch -db sra -query '(poxvirus OR vaccinia OR myxoma OR monkeypox OR mpox) AND \
  ("RNA-seq" OR transcriptome)' | efetch -format runinfo
```

Search terms covered: poxvirus / vaccinia / myxoma / monkeypox / mpox /
orthopoxvirus / cowpox / ectromelia / orf virus, each crossed with RNA-seq,
transcriptome, ribosome profiling, Ribo-seq, polysome, and host shutoff.

## Inclusion Criteria (all required)

1. Public accession available (GEO/SRA/ENA/ArrayExpress).
2. Measures the response for the host for a transcriptomic (or proteomic) response to poxvirus
  Infected and control conditions are both identifiable from metadata.
3. At least 2 biological replicates per condition when they are being used for differential
  expression.
4. Raw counts, FASTQ-derivable counts, or also usable processed result tables have to exist.
5. Metadata that is sufficient to define a clean contrast.

## Exclusion Criteria (any disqualifies from primary meta-analysis)

- No control condition.
- No host expression data (viral-genome-only sequencing).
- Only viral transcript quantification.
- Insufficient metadata to define infection status.
- Single replicate per condition.
- Non-comparable assay (e.g., specialized 5'-leader sequencing).
- Non-mappable host species for translation-factor orthologs.

## Data Tiering


| Tier  | Definition                             | Use                                      |
| ----- | -------------------------------------- | ---------------------------------------- |
| **A** | Raw/processed counts allow our own DGE | DGE, meta-analysis, ML, modules          |
| **B** | Only published DE result tables        | Directional validation, rank correlation |
| **C** | Cannot be cleanly harmonized           | Contextual / hypothesis support only     |


## Outcome

- **Tier A (4 studies):** GSE278320, GSE284044, GSE287860, GSE288000. There were all entered in the primary random-effects meta-analysis (10 harmonized contrasts).
- **Tier B (1 study):** GSE185520. This is external directional validation only.
- **Tier C (excluded from quantitative integration):** GSE288433, GSE329296, GSE228963, GSE137757, GSE280044, GSE228918, GSE259380. This is retained in the registry for transparency and future expansion. See `docs/excluded_datasets.csv`.

This separation is more of a deliberate one and it forces messy or even forces there to be a non-comparable public data into a single model that would then cause an inflation in the apparent significance. Honest tiering is one of the core design principles that are going to be implemented in this atlas.