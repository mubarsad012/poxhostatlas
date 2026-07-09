# Progress log

## 2026-06-24

Today, I completed the initialization of the repository and I wrote a draft of the introduction/background section that is due on Friday. Tomorrow, I also have to make a final image of what I think the research question is going to be so that I can also have this submitted in a timely manner and recieve any necessary feedback. But as of now I am thinking for the central, genome-wide question to be: which host genes and pathways are reproducibly remodeled (what I mean by this is like a consisted reshaping by the virus) across independent poxvirus trascriptome datasets; furthermore, which host hubs are targetted the most? I will decide the exact coherent question tomorrow. 

## 2026-06-25

Today, I have written the research question that I want this project to work toward, key objectives of the project, and made some hypotheses. Furthermore, I have laid out the dataset selection (registry, criteria for exclusion/inclusion, tiering, and expansion targets). 

## 2026-06-27

Today, I built the data-acquisition and also I built the preprocessing pipeline then revieweing the code and making fixes that was making it fail to run with the help of approved Claude AI. Also this pipeline was made for the vaccina anchor dataset GSE278320 (GEO fetach, count harmonization, gene annotation). Also today I also did the start of the [analysis.py](http://analysis.py) which is going to be the reproducible driver. 

## 2026-06-29

Today, the main thing that I worked on was that I wrote and now we also just ran the differential expression with the PyDESeq2 (VacV vs mock) and also then I extracted the translation-factor/RNA-helicase view. 

## 2026-07-01

Today what I did was that I mainly was working on and completed the generation of the primary figures (volcano plot and the heatmap) and also today I worked on adding a strict-title sensitivity analysis for the purpose of robustness to GEO metadata quirks.

## 2026-07-03

Today what I did was that I added a host-only mechanistic curation layer to my project (this was GENCODE biotype filtering and curated host targets) and also added a STRING-based prioritization.

## 2026-07-05

Today I brought in some of the additional datasets. For example what I did was that I added a cross-dataset pilot and expanded differential expression across the Myxoma effector studies and also the a vaccinia time course.

## 2026-07-07

Today what I did was that I was working on figuringi out how to run the genome-wide random-effects meta-analysis (DerSimonian-Laird) and this was for across the independent studies.

## 2026-07-08

Today I was working on adding a leave-one-study-out robustness and also to add a heterogeneity classification to the meta-analysis that I currently have which I think should make this analysis stronger. 

## 2026-07-09

Today I was working on building the advanced and integrative synthesis (in this specifically working on the host-background robustness, temporal dynamics, cross-study rank correlation, and also the composite evidence).