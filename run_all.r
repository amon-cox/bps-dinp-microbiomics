# R script to run all analyses and figure generations in sequence.

## core packages
library(tibble)
library(readr)
library(dplyr)
library(tidyr)
library(purrr)

## visualizations
library(ggplot2)
library(cowplot)

## stats
library(caret)
library(vegan)

## prepare the metadata ahead of running scripts
metadata <- read_csv(file.path("data", "sample_metadata.csv")) |>
    mutate(
        age_week = factor(age_week, levels = c("wk07", "wk12", "wk17", "wk22", "wk40"))
    ) # exposure colors: blue / inoculum, green / control, reddish-orange / bps, reddish-pink / dinp

## call analytical scripts in sequence
#source(file.path("scripts", "1_preprocess.r"))
