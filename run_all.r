# R script to run all analyses in and figure generations in sequence.
library(tibble)
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(caret)
library(ggplot2)
library(cowplot)

metadata <- read_csv(file.path("data", "sample_metadata.csv")) |>
    mutate(
        exposure = factor(exposure, levels = c("inoculum", "control", "bps", "dinp")),
        age_week = factor(age_week, levels = c("wk07", "wk12", "wk17", "wk22", "wk40"))
    )

# exposure colors: blue for inoculum, green for control, reddish-orange for bps, reddish-pink for dinp

source(file.path("scripts", "1_preprocess.r"))
