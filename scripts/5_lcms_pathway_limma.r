# 5) LC-MS putative pathways analyzed by limma
if(!exists("lcms_pathways", mode = "list")) {
    lcms_pathways <- list.files(path = file.path("data", "processed", "lcms"),
        pattern = "^4_",
        full.names = TRUE
    ) |>
    set_names(~ basename(.x) %>% sub("^4_(.*)\\.tsv$", "\\1", .)) |> 
    map(.f = read_tsv)
}

## limma workflow by age-specific dataset
### establish function for applying limma
limma_by_age <- function(df, df_name) {

    #### subset metadata
    md <- metadata |>
        dplyr::filter(sample %in% df$sample, culture_day != "d0") |>
        droplevels()
    
    #### specify storage object
    results_list <- list()

    #### separate data for day-by-day comparisons
    for (day in levels(md$culture_day)) {

        md_day <- md |>
            dplyr::filter(culture_day == day) |>
            droplevels() |>
            mutate(exposure = relevel(exposure, ref = "control"))

        df_day <- dplyr::filter(df, df$sample %in% md_day$sample) |>
            column_to_rownames("sample") |>
            t() |>
            as.data.frame()
        
        df_day <- df_day[, match(md_day$sample, colnames(df_day)), drop = FALSE]

        ##### create design matrix for current time point
        design_matrix <- model.matrix(~ exposure, data = md_day)

        ##### prepare correlations and contrasts
        correlations <- duplicateCorrelation(
            df_day,
            design = design_matrix,
            block = md_day$replicate
        )

        ##### fit the linear models
        fit <- lmFit(
            object = df_day,
            design = design_matrix,
            block = md_day$replicate,
            correlation = correlations$consensus.correlation
        ) |>
            eBayes()

        ##### extract the results
        results <- topTable(
            fit = fit,
            coef = paste0("exposure", setdiff(levels(md_day$exposure), "control")),
            adjust = "BH",
            number = Inf,
            sort.by = "p"
        ) |>
            rownames_to_column("pathway")

        results_list[[day]] <- results
    }

    #### combine results from each day then export
    limma_output <- bind_rows(results_list, .id = "culture_day")

    write_tsv(limma_output, file.path("output", "tables", paste0("5_lcms_pathways_limma_", df_name, ".tsv")))

    sig_limma_res <- limma_output |>
        dplyr::filter(adj.P.Val < 0.05) |>
        as_tibble()
    
    return(sig_limma_res)
}

### run limma function on data list
set.seed(123)

sig_limma_results <- imap(
    .x = lcms_pathways,
    .f = limma_by_age
)

## summary tables of overlaps between ages
sig_summaries <- imap(
    .x = sig_limma_results,
    .f = \(df, i) df |>
        mutate(culture_day = "total_unique") |>
        bind_rows(df) |>
        select(culture_day, pathway) |>
        distinct() |>
        group_by(culture_day) |>
        summarize(count = n(), .groups = "drop")
) |>
    bind_rows(.id = "dataset") |>
    pivot_wider(names_from = "culture_day", values_from = "count", values_fill = 0)

write_csv(
    sig_summaries,
    file.path("output", "tables", "5_lcms_limma_sig_pathway_counts.csv")
)

## summarize overlapping features across ages per cohort
### establish function for overlap comparisons
overlapping_paths <- function(ls) {
    path_list <- map(ls, .f = \(x) pull(x, pathway) |> unique())
    sapply(path_list, function(x) sapply(path_list, function(y) length(intersect(x, y))))
}

### BPS high
sig_limma_results %>%
    purrr::keep(startsWith(names(.), "BPS_high_") & !endsWith(names(.), "wk40")) |>
    overlapping_paths() |>
    write.csv(file.path("output", "tables", "5_lcms_limma_sig_pathway_overalps_BPS_high.csv"))

### DINP high
sig_limma_results %>%
    purrr::keep(startsWith(names(.), "DINP_high_") & !endsWith(names(.), "wk40")) |>
    overlapping_paths() |>
    write.csv(file.path("output", "tables", "5_lcms_limma_sig_pathway_overalps_DINP_high.csv"))

### BPS low
sig_limma_results %>%
    purrr::keep(startsWith(names(.), "BPS_low_")) |>
    overlapping_paths() |>
    write.csv(file.path("output", "tables", "5_lcms_limma_sig_pathway_overalps_BPS_low.csv"))

### DINP low
sig_limma_results %>%
    purrr::keep(startsWith(names(.), "DINP_low_")) |>
    overlapping_paths() |>
    write.csv(file.path("output", "tables", "5_lcms_limma_sig_pathway_overalps_DINP_low.csv"))
