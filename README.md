# PoxHostAtlas — a reproducible cross-study analysis of host gene-expression remodeling during poxvirus infection (work in progress)

PoxHostAtlas, as of my current plan of action, is going to be a computational secondary-analysis project that asks: **genome-wide and unbiased, which host genes and pathways are reproducibly remodeled the most across independent public poxvirus transcriptome datasets**. And to start this off, I will begin with the prototypical poxvirus, which is the vaccinia virus.

> **Current Status Check: early / in progress.** I am planning to grow this README and it will grow as the project develops. The aim, the background, and also the plan will be located in `docs/`.

## Aim (mentor guidance also considered)

The aim of this project will be to gather the independent public poxvirus host-response RNA-seq. After that I hope to then harmonize differential expression and then identify the host genes or pathways that are remodeled reproducibly across labs and contexts. Furthermore, I will also check the host hubs that are targeted the most consistently. At this point we are not deliberately pre-narrowing onto one protein family as per the guidance of my mentor. RNA-helicase/translation remodeling will be examined as one sub-question within the broader map that we have. 

- Background & motivation: `[docs/01_INTRODUCTION.md](docs/01_INTRODUCTION.md)`
- Research question & plan: `[docs/02_RESEARCH_QUESTION.md](docs/02_RESEARCH_QUESTION.md)`
- Data sources: `[docs/SOURCES.md](docs/SOURCES.md)`

## Setup

```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

## Progress log

Just to remain organized and track everything, the progress log for the current whereabouts of the project will be in `docs/progress/`.

## Claude AI

Claude Code (Anthropic) will be utilized within this project as a supplementary tool for code verification and formatting during the preparation of this research project and also the accompanying figures. Approval was recieved from my mentor, as OpenBio SRI has directed students to do. 