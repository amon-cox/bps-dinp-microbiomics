# bps-dinp-microbiomics

This repository contains the data, analytical code, and results for a study on age‑dependent microbiome and metabolomic responses to BPS and DINP in murine fecal cultures.

Working title: Investigating the Influence of Host Age on Murine Microbial Responses to Dietary Contaminants.

All analyses and code in this repository were produced by the authors for an original research manuscript that has not yet been peer‑reviewed. The code contained here is a reorganization of the original analytical code, incorporating updates for reproducibility and compatability between operating systems. Please do not cite this repository as the definitive source; a DOI and manuscript reference will be added after peer review. Contents may change until the manuscript is submitted, reviewed, and approved for publishing, whereupon this repository will be updated with the final citation and archived with Zenodo.

## Authorship and affiliation

Amon Cox\*, Evelyn Callaway, and Arul Jayaraman\*\*

Artie McFerrin Department of Chemical Engineering, Texas A&M University, College Station, Texas, United States of America, 77843

\*Repository correspendence: <amoncox@gmail.com> | GitHub: amon-cox

\*\*Research correspondence: <arulj@tamu.edu>

## License

- Code in the `scripts/` directory is licensed under the MIT License (see `LICENSE`).
- Data, figures, and objects in the `data/` and `output/` directories and `Supplementary Material.pdf` are licensed under the Creative Commons Attribution 4.0 International license (CC BY 4.0; see `LICENSE-CC-BY-4.0`).

## Abstract

Host age can influence microbiome structure and function, but its impact on microbial responses to dietary contaminants remains unclear. Fecal microbiota from female C57BL/6J mice (7-40 weeks) were cultured in vitro with bisphenol S or diisononyl phthalate at 10 or 100 μM and profiled by 16S rRNA gene sequencing and untargeted LC-MS metabolomics. We hypothesized that communities derived from younger mice would show the greatest dose-dependent disruption. In contrast, the metabolomic profiles from mature (22-week) and aging (40-week) donors exhibited the most pronounced, age-dependent alterations following 100 μM BPS or DINP, whereas 10 μM exposures produced minimal changes beyond culture-time effects. Despite robust metabolic shifts, no significant differences in community composition or individual taxa were detected, indicating functional plasticity without dysbiosis. Inspecting metabolomic and predicted functional profiles revealed age- and exposure-responsive metabolic pathways, particularly short-chain fatty acid and amino acid metabolism for DINP, and indole alkaloid and lipid-related metabolism for BPS, some of which overlap with known bisphenol and phthalate effects. This study suggests that host age affects microbiome metabolic susceptibility to dietary plastics and identifies a small set of metabolic endpoints that may serve as sensitive markers of microbiome vulnerability in future studies using human-relevant exposure levels.

## Repository Structure

- `data/`
  - `original/` – input data and sample metadata (read-only).
    - `16s/` – relative abundances and taxanomic metadata from 16S rRNA sequencing vendor.
    - `lcms/` – tabulated LC-MS peak intensities from LC-MS platform.
  - `processed/` – filtered, quality-controled, reformatted data used as inputs to analysis scripts.
    - `16s/` – ASV tables and pathway abundances predicted by Tax4Fun2 on MicrobiomeAnalyst.
    - `lcms/` – normalized LC‑MS feature tables and putative pathway annotations.

- `output/`
  - `16s_Maaslin/` – MaAsLin outputs for 16S analyses (model summaries, coefficients).
  - `figures/` – rendered figures used in the manuscript and Supplementary Material.
  - `tables/` – intermediate and final results tables (e.g., PCA, limma, PERMANOVA).

- `scripts/` – analysis scripts, numerated by run order.

- `renv/` – project‑local package library and metadata for reproducible R environments (automatic).

- `run_all.R` – convenience script that runs the full analysis pipeline
  to reproduce key tables and figures from the manuscript.

- `supplementary.qmd` – Quarto source for the Supplementary Material.
- `Supplementary Material.pdf` – rendered Supplementary Material corresponding to the manuscript.

## How to run

To reproduce the contents of `data/processed/` and `output/`, open an R session from the project root directory and run the following:

```r
install.packages("renv")
renv::restore()
source("run_all.r")
```
