# 1) Prepare all datasets for analysis.
## organize LC-MS datasets
### call in LC-MS datasets as list
lcms_peaks <- list.files(path = file.path("data", "original", "lcms"),
        pattern = "^lcms_peaks_",
        full.names = TRUE
    ) |>
    set_names(~ basename(.x) %>% sub("\\.csv$", "", .)) |> 
    map(.f = read_csv)

### remove duplicated rows (artefact) and transform to log2
lcms_peaks_log2 <- map(
    .x = lcms_peaks,
    .f = function(.df) distinct(.df) |> mutate(across(where(is.numeric), log2))
)

### reorganize data by compound ($cohort), dose level, and age
lcms_peaks_long <- lcms_peaks_log2 |>
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
    filter(!is.na(intensity)) |>
    group_by(cohort, dose, age_week)

lcms_key <- group_keys(lcms_peaks_long) |>
    mutate(group_id = paste(cohort, dose, age_week, sep = "_"))

lcms_peaks_grouped <- lcms_peaks_long |>
    group_map(~ .x |>
            ungroup() |>
            pivot_wider(
                names_from = mz__rtMin,
                values_from = intensity
            )
    )

names(lcms_peaks_grouped) <- pull(lcms_key, group_id)

## organize 16S datasets
### call in 16S datasets
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
    filter(!is.na(relative_abundance)) |>
    group_by(cohort, dose, age_week)

seq_key <- group_keys(seq_ab_long) |>
    mutate(group_id = paste(cohort, dose, age_week, sep = "_"))

### remove low-variance features by group then export
seq_noLowVar_grouped <- seq_ab_long |>
    group_map(~ {
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
            rownames_to_column("sample") |>
            as_tibble()
    })

names(seq_noLowVar_grouped) <- pull(seq_key, group_id)

## condcut PCA to gauge variances and identify outliers
### prepare workflow for PCA, outlier detection, and export
outliers_full <- list()

pca_and_outliers <- function(.df, .df_name, .ome) {

    #### conduct PCA and calculate variance
    pca <- .df |>
        column_to_rownames("sample") |>
        prcomp(center = TRUE, scale. = TRUE)

    variance <- pca$sdev^2 / sum(pca$sdev^2)

    #### export PCA scores
    pca$x |>
        as.data.frame() |>
        rownames_to_column("sample") |>
        write_csv(file.path("output", "tables", paste0("1_pca_scores_", .df_name, ".csv")))

    #### determine z-scores and outliers based off of PCs 1 & 2
    pc_z_scores <- scale(pca$x[, 1:2])

    outliers <- row.names(pca$x)[abs(pc_z_scores[, 1]) > 3 | abs(pc_z_scores[, 2]) > 3]

    #### remove outliers from current dataset then export
    if (length(outliers) != 0) {
        outliers_full[[paste0(.ome, "_", .df_name)]] <<- outliers
    }
    
    data_no_outliers <- filter(.df, !sample %in% outliers)

    write_csv(
        x = data_no_outliers,
        file = file.path("data", "processed", .ome, paste0("1_", .df_name, ".csv"))
    )

    invisible(NULL)
}

### run the PCA and outlier function on the datasets
#### LC-MS data
map2(
    .x = lcms_peaks_grouped,
    .y = names(lcms_peaks_grouped),
    .f = \(df, nm) pca_and_outliers(.df = df, .df_name = nm, .ome = "lcms")
)

#### 16S data
map2(
    .x = seq_noLowVar_grouped,
    .y = names(seq_noLowVar_grouped),
    .f = \(df, nm) pca_and_outliers(.df = df, .df_name = nm, .ome = "16s")
)

### format outliers_full to accomodate data type (16s vs lcms) and dataset (cohort, dose, age_week)
outliers_full |>
    bind_rows(.id = "dataset") |>
    pivot_longer(cols = everything(), names_to = "dataset", values_to = "outlier") |>
    separate_wider_regex(
        dataset,
        patterns = c(
            data_type = "^[^_]+",
            "_",
            dataset_name = ".*"
        )
    ) |>
    distinct() |>
    left_join(
        y = select(metadata, sample, cohort, dose, age_week),
        by = c("outlier" = "sample")
    ) |>
    arrange(data_type, dose, cohort, age_week) |>
    relocate(data_type, outlier, dose, cohort, age_week, dataset_name) |>
    write_csv(file.path("output", "tables", "1_outliers.csv"))
