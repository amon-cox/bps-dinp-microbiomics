# 11) MaAslin2 comparison of compositional baseline for mice aged 22 and 40 weeks, using DINP cohort 16S data.
if(!exists("rel_ab_processed", mode = "list")) {
    rel_ab_processed <- list.files(path = file.path("data", "processed", "16s"),
        pattern = "^1_",
        full.names = TRUE
    ) |>
        set_names(~ basename(.x) %>% sub("^1_(.*)\\.csv$", "\\1", .)) |> 
        map(.f = read_csv)
}

## subset data to inocula (day 0) samples and average the week 22 samples by mouse
rel_ab_d0 <- rel_ab_processed |>
    bind_rows() |>
    left_join(y = metadata, by = "sample") |>
    dplyr::filter(culture_day == "d0") |>
    select(sample, age_week, replicate, starts_with("ASV")) |>
    summarize(across(starts_with("ASV"), mean), .by = c("age_week", "replicate")) |>
    mutate(subject = paste(age_week, replicate, sep = "_"))

## run MaAsLin2 comparing inocula between ages
Maaslin2(
    input_data = rel_ab_d0 |> select(subject, starts_with("ASV")) |> column_to_rownames("subject"),
    input_metadata = column_to_rownames(rel_ab_d0, "subject"),
    min_abundance = 0.0001, # setting minimum abundance (0.01%)
    min_prevalence = 0.1, # must be prevalent in 10% of samples for consideration
    min_variance = 0,
    normalization = "NONE", # total-sum scaling
    transform = "AST", # recommended arcsine square root for proportional data
    fixed_effects = "age_week",
    random_effects = "replicate",
    analysis_method = "LM", # recommended "LM" method for percent abundance data
    reference = "age_week,wk22",
    correction = "BH", # Bonferroni correction for q values
    output = file.path(
        "output", 
        "16s_Maaslin", 
        "11_16s_Maaslin2_inocula_wk22_wk40"
    )
)
