# 1) Prepare all datasets for analysis.
## process LC-MS datasets
### call in LC-MS datasets as list
lcms_peaks <- list.files(path = file.path("data", "original", "lcms"),
        pattern = "^lcms_peaks_",
        full.names = TRUE
    ) |>
    set_names(~ basename(.x) %>% sub("\\.csv$", "", .)) |> 
    map(.f = read_csv)

### remove duplicated rows (artefact) and transform to log2
lcms_peaks_processed <- map(
    .x = lcms_peaks,
    .f = function(.df) distinct(.df) |> mutate(across(where(is.numeric), log2))
)

### reorganize data by compound ($cohort), dose level, and age
lcms_peaks_long <- lcms_peaks_processed |>
    bind_rows() |>
    pivot_longer(
        cols = -mz__rtMin,
        names_to = "sample",
        values_to = "intensity"
    ) |>
    left_join(
        y = select(metadata, sample, cohort, dose, age_week),
        by = "sample"
    ) |>
    filter(!is.na(intensity))

### export by the groups for distinct, age-based datasets
lcms_peaks_long |>
    group_by(cohort, dose, age_week) |>
    group_walk(~ .x |>
            ungroup() |>
            pivot_wider(
                names_from = mz__rtMin,
                values_from = intensity
            ) |>
            write_csv(
                file.path(
                "data", "processed", "lcms",
                paste0(.y$cohort, "_", .y$dose, "_", .y$age_week, ".csv")
                )
            )
    )

## process 16S dataset
### call in 16S dataset
seq_rel_ab <- read_csv(file.path("data", "original", "16s", "sequencing_relative_abundances.csv"))

### split data by groups
seq_ab_long <- seq_rel_ab |>
    pivot_longer(
        cols = -ASV,
        names_to = "sample",
        values_to = "relative_abundance"
    ) |>
    left_join(
        y = select(metadata, sample, cohort, dose, age_week),
        by = "sample"
    ) |>
    filter(!is.na(relative_abundance))

### remove low-variance features by group then export
seq_ab_long |>
    group_by(cohort, dose, age_week) |>
    group_walk(~ {
        wide <- .x |>
            ungroup() |>
            pivot_wider(
                names_from = ASV,
                values_from = relative_abundance
            )
        
        nz <- nearZeroVar(
            wide |> select(-sample),
            saveMetrics = FALSE
        )

        wide_noZeroVar <- wide |> 
            column_to_rownames("sample") |>
            select(-all_of(nz)) |>
            rownames_to_column("sample")

        write_csv(
            x = wide_noZeroVar,
            file = file.path(
                "data", "processed", "16s",
                paste0(.y$cohort, "_", .y$dose, "_", .y$age_week, ".csv")
            )
        )
    })
