# 10) Differential analysis on predicted pathways from 16S data
if(!exists("t4f2_path_ls", mode = "list")) {
    t4f2_path_ls <- list.files(path = file.path("data", "processed", "16s"),
        pattern = "^6_(.*)\\_pathways.csv",
        full.names = TRUE
    ) |>
    set_names(~ basename(.x) %>% sub("^6_(.*)\\_pathways.csv$", "\\1", .)) |> 
    map(.f = read_csv)
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

    write_tsv(limma_output, file.path("output", "tables", paste0("10_16s_pathways_limma_", df_name, ".tsv")))

    sig_limma_res <- limma_output |>
        dplyr::filter(adj.P.Val < 0.05) |>
        as_tibble()
    
    return(sig_limma_res)
}

### run limma function on data list
set.seed(123)

sig_limma_16s_results <- imap(
    .x = t4f2_path_ls,
    .f = limma_by_age
)
