# R script to run all analyses and figure generations in sequence.
## script to load packages
source(file.path("scripts", "0_setup.r"))

## call analytical scripts in sequence
source(file.path("scripts", "1_preprocess.r"))
source(file.path("scripts", "2_lcms_profiling.r"))
source(file.path("scripts", "3_lcms_limma.r"))
source(file.path("scripts", "4_lcms_pathway_annotation.r"))
source(file.path("scripts", "5_lcms_pathway_limma.r"))
source(file.path("scripts", "6_pathway_comparison."))
source(file.path("scripts", "7_16s_composition.r"))
source(file.path("scripts", "8_16s_diversity.r"))
source(file.path("scripts", "9_16s_maaslin.r"))
source(file.path("scripts", "10_16s_pathway_limma.r"))
source(file.path("scripts", "11_16s_compare_inocula."))
source(file.path("scripts", "12_lcms_compare_inocula.r"))
source(file.path("scripts", "13_age_and_exposure.r"))

## produce supplementary file
quarto::quarto_render(
    input = "supplementary.qmd",
    output_format = "pdf",
    output_file = "Supplementary Material.pdf"
)
