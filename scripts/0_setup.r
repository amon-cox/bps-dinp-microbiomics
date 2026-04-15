# 0) load all necessary packages and sample metadata
## data management
library(tibble)
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(forcats)

## visualizations
library(ggplot2)
library(cowplot)

## stats
library(caret)
library(vegan)
library(limma)
library(FSA)
library(Maaslin2)
library(pairwiseAdonis) # remotes::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

## prepare the metadata ahead of running scripts
metadata <- read_csv(file.path("data", "sample_metadata.csv")) |>
    mutate(
        age_week = factor(age_week, levels = c("wk07", "wk12", "wk17", "wk22", "wk40")),
        exposure = factor(exposure, levels = c("inoculum", "control", "bps", "dinp"))
    )

## prepare aesthetics for plotting
exposure_color_palette <- viridisLite::turbo(length(levels(metadata$exposure)), begin = 0.1, end = 0.9)
names(exposure_color_palette) <- levels(metadata$exposure)
