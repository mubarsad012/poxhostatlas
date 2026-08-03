# Sources

## Dataset

- GEO accession: GSE278320
- GEO accession page: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278320>
- GEO supplementary file list: <https://ftp.ncbi.nlm.nih.gov/geo/series/GSE278nnn/GSE278320/suppl/filelist.txt>
- GEO family SOFT file: <https://ftp.ncbi.nlm.nih.gov/geo/series/GSE278nnn/GSE278320/soft/GSE278320_family.soft.gz>

## Primary Publication

- Park C, Ferrell AJ, Meade N, Shen PS, Walsh D, and colleagues. Distinct non-canonical translation initiation modes arise for specific host and viral mRNAs during poxvirus-induced shutoff. Nature Microbiology. 2025;10:1535-1549. doi:10.1038/s41564-025-02009-4.
- Article page: <https://www.nature.com/articles/s41564-025-02009-4>
- PubMed Central record: <https://pmc.ncbi.nlm.nih.gov/articles/PMC12305801/>

## Analysis Software

- PyDESeq2 package: <https://pypi.org/project/pydeseq2/>
- PyDESeq2 workflow documentation: <https://pydeseq2.readthedocs.io/en/latest/auto_examples/plot_minimal_pydeseq2_pipeline.html>
- DESeq2 statistical method: Love MI, Huber W, Anders S. Genome Biology. 2014. doi:10.1186/s13059-014-0550-8.
- GENCODE v44 annotation directory: <https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/>
- GENCODE v44 GTF used for host biotype curation: <https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz>
- STRING API documentation: <https://string-db.org/help/api/>
- STRING 2023 database reference: Szklarczyk et al. Nucleic Acids Research. 2023. doi:10.1093/nar/gkac1000.

## Poxvirus Transcriptome-Over-Time Literature (mentor-flagged, 2026-06)

These four papers were flagged by the project mentor as the key "poxvirus
transcriptome over time" references. Two profile the **viral** transcriptome
(temporal/kinetics context); two profile the **host** response (expansion
targets). See `docs/SUPERVISOR_UPDATE_2026-06-21.md` for the usability triage.

- Tombácz et al. Time-Course Transcriptome Profiling of a Poxvirus Using Long-Read Full-Length Assay. Pathogens. 2021;10(8):919. PMID 34451383. <https://pmc.ncbi.nlm.nih.gov/articles/PMC8398953/> — long-read **viral** transcriptome kinetics (context only). Data: ENA **PRJEB26434** (Characterization of the Vaccinia virus transcriptome) and **PRJEB26430** (Dynamic characterization of the Vaccinia virus transcriptome) — long-read subreads/direct-RNA, not host count matrices.
- Tombácz et al. Dynamic transcriptome profiling dataset of vaccinia virus obtained from long-read sequencing techniques. GigaScience. 2018;7(12):giy139. PMID 30476066. ENA **PRJEB26434**. <https://pmc.ncbi.nlm.nih.gov/articles/PMC6290886/> — long-read **viral** transcriptome data note (context only).
- Yang Z, Cao S, Martens CA, Porcella SF, Xie Z, Ma M, Shen B, Moss B. Deciphering Poxvirus Gene Expression by RNA Sequencing and Ribosome Profiling. J Virol. 2015;89(13):6874-6886. PMID 25903347; doi:10.1128/JVI.00528-15. <https://pmc.ncbi.nlm.nih.gov/articles/PMC4468498/> — bulk **host** RNA-seq + ribosome profiling time course. Data: BioProject **PRJNA280609** / SRA **SRP056975** (raw FASTQ only; **no GEO processed matrix**) — 3 RNA-seq + 12 ribosome-profiling runs. A ready-to-run download/align/quantify pipeline is vendored at `data/external/yang_moss_2015/` (needs network + tools). High-priority host expansion target; mock-vs-infected contrast pending verification of sample titles.
- Matía A, McCarthy F, Woosley H, et al. Spatio-temporal analysis of Vaccinia virus infection and host response dynamics using single-cell transcriptomics and proteomics. bioRxiv. 2024. doi:10.1101/2024.01.13.575413. <https://www.biorxiv.org/content/10.1101/2024.01.13.575413v1> — **host** single-cell RNA-seq + proteomics. **Integrated as external validation** (`scripts/24_matia_external_validation.py`): supplementary single-cell metadata (HeLa + BJ5ta) and HeLa proteomics are vendored under `data/external/matia2024/`. The raw scRNA-seq count matrices' accession is not currently indexed in public GEO/SRA searches (preprint; data-availability not resolvable from the deposited text), so validation uses the supplementary tables directly. Preprint, not peer-reviewed.

## Mechanistic Literature

- DHX29/eIF3 scanning mechanism: Pisareva VP, Pisarev AV. DHX29 and eIF3 cooperate in ribosomal scanning on structured mRNAs during translation initiation. RNA. 2016. doi:10.1261/rna.057851.116.
- Myxoma/DHX helicase screen: Identification of host DEAD-box RNA helicases that regulate cellular tropism of oncolytic Myxoma virus in human cancer cells. Scientific Reports. 2017. doi:10.1038/s41598-017-15941-1.

## Primary Sample Selection

The primary analysis compares parental HAP1 total-fraction count files.
Supplementary count filenames are used as the fraction source of truth, and
GEO title/file disagreements are retained in `data/processed/sample_manifest.csv`.
