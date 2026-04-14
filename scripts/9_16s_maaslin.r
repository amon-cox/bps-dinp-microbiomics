# 9) identifying differentially abundant ASVs with MaAsLin2
if(!exists("rel_ab_processed", mode = "list")) {
    rel_ab_processed <- list.files(path = file.path("data", "processed", "16s"),
        pattern = "^1_",
        full.names = TRUE
    ) |>
        set_names(~ basename(.x) %>% sub("^1_(.*)\\.csv$", "\\1", .)) |> 
        map(.f = read_csv)
}

## perform MaAsLin2 for each cohort x dose x age
imap(
    .x = rel_ab_processed,
    .f = \(df, nm) {
        ### 
        df_no_inocula <- df |>
            left_join(select(metadata, sample, exposure), by = "sample") |>
            dplyr::filter(exposure != "inoculum") |>
            select(-exposure)

        ### subset metadata
        md_subset <- metadata |>
            dplyr::filter(sample %in% df_no_inocula$sample) |>
            droplevels()
        
        ### run Maaslin2 function
        Maaslin2(
            input_data = column_to_rownames(df, "sample"),
            input_metadata = column_to_rownames(md_subset, "sample"),
            min_abundance = 0.0001, # setting minimum abundance (0.01%)
            min_prevalence = 0.1, # must be prevalent in 10% of samples for consideration
            min_variance = 0,
            normalization = "NONE", # total-sum scaling
            transform = "AST", # recommended arcsine square root for proportional data
            fixed_effects = c("exposure", "culture_day", "exposure:culture_day"),
            random_effects = "replicate",
            analysis_method = "LM", # recommended "LM" method for percent abundance data
            reference = c("exposure,control", "culture_day,d2"),
            correction = "BH", # Bonferroni correction for q values
            output = file.path("output", "16s_Maaslin", paste0("9_16s_Maaslin2_", nm))
        )
    }
)
