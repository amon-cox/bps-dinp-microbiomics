# 2) PERMANOVA and Principle Coordinates Analysis for LC-MS data
## call in processed data
lcms_peaks_processed <- list.files(path = file.path("data", "processed", "lcms"),
        pattern = "^1_",
        full.names = TRUE
    ) |>
    set_names(~ basename(.x) %>% sub("^1_(.*)\\.csv$", "\\1", .)) |> 
    map(.f = read_csv)

## establish distance matrices
set.seed(123)

lcms_distances <- map(
    .x = lcms_peaks_processed,
    .f = \(x) column_to_rownames(x, "sample") |> 
        scale(center = TRUE, scale = TRUE) |> 
        vegdist(method = "euclidean")
)

## function for performing PERMANOVA
permanova_by_set <- function(ds, ds_name) {

    ### subset metadata
    md <- dplyr::filter(metadata, sample %in% row.names(as.matrix(ds)))

    ### global PERMANOVA
    permanova_global <- adonis2(
        ds ~ exposure + culture_day,
        data = md,
        strata = md$replicate
    )

    ### PERMANOVA by day
    days <- setdiff(unique(md$culture_day), "d0")

    permanova_daily <- lapply(days, function(day_i) {
        idx <- md$culture_day == day_i
        samples_day <- md$sample[idx]
        replicates_day <- md$replicate[idx]

        ds_m <- as.matrix(ds)[samples_day, samples_day]
        ds_day <- as.dist(ds_m)
        md_day <- md[idx, , drop = FALSE]

        res <- adonis2(
            ds_day ~ exposure,
            data = md_day,
            strata = replicates_day
        )
        rownames(res) <- c("exposure", "Residual", "Total")
        res
    })

    names(permanova_daily) <- days

    global_rows <- tibble(
        PERMANOVA = c("global: exposure"),
        R2 = round(permanova_global$R2[1], 4),
        p_val = round(permanova_global[1, "Pr(>F)"], 4)
    )

    daily_rows <- imap_dfr(
        permanova_daily,
        ~ tibble(
            PERMANOVA = paste0(.y, ": exposure"),
            R2 = round(.x[1, "R2"], 4),
            p_val = round(.x[1, "Pr(>F)"], 4)
        )
    )

    permanova_table <- bind_rows(global_rows, daily_rows)

    write_csv(
        permanova_table,
        file.path("output", "tables", paste0("2_lcms_PERMANOVA_", ds_name, ".csv"))
    )
}

map2(
    .x = lcms_distances,
    .y = names(lcms_distances),
    .f = \(x, y) permanova_by_set(ds = x, ds_name = y)
)
