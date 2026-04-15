# bps-dinp-microbiomics

Code and data for a study on age‑dependent microbiome and metabolomic responses to BPS and DINP in murine fecal cultures. All analyses and code in this repository were developed by the author(s) for an original research manuscript that has not yet been peer‑reviewed.

Working title: Investigating the Influence of Host Age on Murine Microbial Responses to Dietary Contaminants.

## Project status

This repository contains analysis code, data, and results for an unpublished research manuscript currently in preparation for submission. The code contained here is a reorganization of the original analytical code, incorporating updates for reproducibility and compatability between operating systems. Please do not cite this repository as the definitive source; a DOI and manuscript reference will be added after peer review. Contents may change until the manuscript is submitted, reviewed, and approved for publishing.

## Authorship and affiliation

Lead author: Amon Cox (<amoncox@gmail.com> | GitHub: amon-cox)

Additional authors and full affiliations will be added following manuscript submission.

## License

- Code in the `scripts/` directory is licensed under the MIT License (see `LICENSE`).
- Data, figures, and objects in the `data/` and `output/` directories are licensed under the Creative Commons Attribution 4.0 International license (CC BY 4.0; see `LICENSE-CC-BY-4.0`).

## Planned updates

- [x] Add original data and project scaffolding
- [ ] Code restructuring and compatability updates (e.g. file.path() insertions,  reproducibility safeguards)
  - [ ] run_all.r
  - [x] 0_setup.r
  - [x] 1_preprocess.r
  - [x] 2_lcms_profiling.r
  - [x] 3_lcms_limma.r
  - [ ] 4_lcms_pathway_annotation.r
  - [ ] 5_lcms_pathway_overview.r
  - [ ] 6_lcms_pathway_limma.r
  - [x] 7_16s_composition.r
  - [x] 8_16s_diversity.r
  - [x] 9_16s_maaslin.r
  - [ ] 10_compare_pathways.r
  - [ ] 11_compare_inocula.r
- [ ] File guide
- [ ] Upon manuscript acceptance, this repository will be tagged, archived with Zenodo, and updated with the final citation.
