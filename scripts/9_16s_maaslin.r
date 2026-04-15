# 9) identifying differentially abundant ASVs with MaAsLin2
if(!exists("rel_ab_processed", mode = "list")) {
    rel_ab_processed <- list.files(path = file.path("data", "processed", "16s"),
        pattern = "^1_",
        full.names = TRUE
    ) |>
        set_names(~ basename(.x) %>% sub("^1_(.*)\\.csv$", "\\1", .)) |> 
        map(.f = read_csv)
}

## prepare wrapper function for MaAsLin2 by day
maaslin_by_day <- function(df, nm) {

    df_metadata <- metadata |>
        dplyr::filter(
            sample %in% df$sample,
            exposure != "inoculum"
        ) |>
        mutate(culture_day = factor(culture_day))
    
    for (day in levels(df_metadata$culture_day)) {
        ### subset the data by day
        df_day <- df |>
            left_join(select(df_metadata, sample, exposure, culture_day, replicate), by = "sample") |>
            dplyr::filter(culture_day == day) |>
            mutate(replicate = factor(replicate)) |>
            droplevels()
        
        features <- df_day |>
            select(-exposure, -culture_day, -replicate) |>
            column_to_rownames("sample")
        
        meta <- df_day |>
            select(sample, exposure, culture_day, replicate) |>
            column_to_rownames("sample")
        
        ### run Maaslin2 function
        Maaslin2(
            input_data = features,
            input_metadata = meta,
            min_abundance = 0.0001, # setting minimum abundance (0.01%)
            min_prevalence = 0.1, # must be prevalent in 10% of samples for consideration
            min_variance = 0,
            normalization = "NONE", # total-sum scaling
            transform = "AST", # recommended arcsine square root for proportional data
            fixed_effects = "exposure",
            random_effects = "replicate",
            analysis_method = "LM", # recommended "LM" method for percent abundance data
            reference = "exposure,control",
            correction = "BH", # Bonferroni correction for q values
            output = file.path(
                "output", 
                "16s_Maaslin", 
                paste0("9_16s_Maaslin2_", nm, "_", day)
            )
        )
    }
}

## perform MaAsLin2 for each cohort x dose x age, split by day
imap(rel_ab_processed, maaslin_by_day)
