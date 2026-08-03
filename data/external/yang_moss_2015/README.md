# Yang & Moss 2015 (SRP056975) — ingestion notes

Deciphering Poxvirus Gene Expression by RNA Sequencing and Ribosome Profiling.
J Virol 2015; PMID 25903347; doi:10.1128/JVI.00528-15.

- BioProject: **PRJNA280609**  | SRA study: **SRP056975**  | BioSample: SAMN03465474
- Host: human HeLa; virus: Vaccinia virus WR.
- 15 runs = **3 RNA-seq** (SRR1959018, SRR1959019, SRR1959020) + **12 ribosome
  profiling / RPF** (`library_strategy = OTHER`: SRR1959021–SRR1959032).
- **Raw FASTQ only.** GEO has no processed count matrix for this study
  (`https://www.ncbi.nlm.nih.gov/gds/?term=SRP056975` returns nothing), so it
  must be downloaded and quantified.

## ⚠️ Before you quantify: map conditions (REQUIRED)

The ENA file report (`filereport_SRP056975.tsv`) does **not** contain sample
titles, so we cannot yet tell which runs are mock vs which infection time point.
Pull the richer report (adds titles), then fill in `sample_conditions.tsv`:

```
https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP056975&result=read_run&fields=run_accession,experiment_accession,experiment_title,sample_title,library_strategy,fastq_ftp&format=tsv&download=true&limit=0
```

Whether this dataset can contribute an **infected-vs-mock host DGE contrast** to
the meta-analysis depends on whether a mock/uninfected RNA-seq sample exists. If
the RNA-seq runs are all infected time points (e.g. 2/4/8 hpi with no mock), the
dataset still contributes (a) temporal host dynamics and (b) **translation
efficiency** (RPF/RNA) — both axes the mentor cares about — but not a clean
mock contrast. Confirm from the titles and the paper before integrating.

## How to run

1. Install tools: `conda env create -f environment.yml && conda activate srp056975`
   (sra-tools, fastp, cutadapt, bowtie2, STAR, samtools, subread).
2. Build references once (human GRCh38 + GENCODE v44 GTF to match the rest of the
   repo, plus the Vaccinia WR genome so viral reads are absorbed and not
   miscounted as host). URLs are in `quantify_SRP056975.sh`.
3. `bash quantify_SRP056975.sh` — downloads FASTQ, trims, aligns, and writes a
   host count matrix in the same format as `data/processed/counts.csv`.
4. Map conditions in `sample_conditions.tsv`, then a standard PyDESeq2 contrast
   (mirroring `scripts/03_run_pydeseq2.py`) yields the host DGE table that can be
   added as a study to `scripts/13_run_meta_analysis.py`.

Compute note: short-read RNA-seq + Ribo-seq alignment needs ~50–100 GB disk and
a multi-core machine; budget a few hours. This is why it is not run inside the
offline repo pipeline.
