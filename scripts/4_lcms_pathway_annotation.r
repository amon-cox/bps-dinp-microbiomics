# 4) Reorganize LC-MS data by pathways based on putative KEGG annotations
if (!exists("features_annotated")) {
    features_annotated <- read_tsv(file.path("data", "processed", "lcms", "putative_feature_annotations.tsv"))
}

if(!exists("lcms_peaks_processed", mode = "list")) {
    lcms_peaks_processed <- list.files(path = file.path("data", "processed", "lcms"),
        pattern = "^1_",
        full.names = TRUE
    ) |>
    set_names(~ basename(.x) %>% sub("^1_(.*)\\.csv$", "\\1", .)) |> 
    map(.f = read_csv)
}

## join pathway info to peaks data, dilute peak intensities across pathways, then summarize pathway expression
lcms_pathways <- imap(
    .x = lcms_peaks_processed,
    .f = \(df, nm) {
        ### reorganize data to merge samples with putative annotations by LC-MS run
        df_metadata <- left_join(
            x = df,
            y = select(metadata, sample, lcms_run),
            by = "sample"
        ) |>
        pivot_longer(
            cols = -c(sample, lcms_run),
            names_to = "peak",
            values_to = "log2_au"
        )

        df_annotated <- left_join(
            x = df_metadata,
            y = features_annotated,
            by = c("lcms_run", "peak" = "mz__rtMin"),
            relationship = "many-to-many"
        ) |>
        drop_na(pathways)

        ### "weigh" each peak intensity by the number of pathways it contributes to, then sum the pathways
        df_pathways <- df_annotated |>
            separate_longer_delim(pathways, delim = ";") |> # expand each peak by its pathway list
            add_count(sample, peak, name = "cases") |> # create a column that tracks the frequency of each peak's appearance
            rowwise() |> # dilute the peak intensities row by row
            mutate(intensity_au = 2^log2_au / cases) |> # convert back to intensities, divide by frequency
            ungroup() |>
            group_by(sample, pathways) |>
            summarize( # sum pathway intensity and convert back to log2
                pathway_intensity_au = log2(sum(intensity_au)),
                .groups = "drop"
            ) |> 
            pivot_wider(
                names_from = pathways,
                values_from = pathway_intensity_au
            )
        
        ### export pathways data
        write_tsv(
            df_pathways,
            file.path("data", "processed", "lcms", paste0("4_", nm, ".tsv"))
        )

        return(df_pathways)
    }
)
