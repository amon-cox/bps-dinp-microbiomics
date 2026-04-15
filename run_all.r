# R script to run all analyses and figure generations in sequence.

## script to load packages
source(file.path("scripts", "0_setup.r"))

## call analytical scripts in sequence
source(file.path("scripts", "1_preprocess.r"))
source(file.path("scripts", "2_lcms_profiling.r"))
source(file.path("scripts", "3_lcms_limma.r"))

source(file.path("scripts", "7_16s_composition.r"))
source(file.path("scripts", "8_16s_diversity.r"))
source(file.path("scripts", "9_16s_maaslin.r"))
