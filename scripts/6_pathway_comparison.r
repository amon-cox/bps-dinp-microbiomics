# 6) CCA on LC-MS log2 peak intensities and 16S predicted pathway counts
if(!exists("lcms_pathways", mode = "list")) {
    lcms_pathways <- list.files(path = file.path("data", "processed", "lcms"),
        pattern = "^4_",
        full.names = TRUE
    ) |>
    set_names(~ basename(.x) %>% sub("^4_(.*)\\.tsv$", "\\1", .)) |> 
    map(.f = read_tsv)
}

t4f2_pathways <- read_csv(file.path("data", "processed", "16s", "MicrobiomeAnalyst", "sequencing_Tax4Fun2_pathways.csv"))

## filter 16S paths data to remove relevant outlier samples.
outliers <- read_csv(file.path("output", "tables", "1_outliers.csv"))
outliers_16s <- dplyr::filter(outliers, data_type == "16s")

t4f2_path_ls <- t4f2_pathways |>
    mutate(dataset = str_match(sample, "^([^-]+)(?:_[^-]+){3}$") %>% .[,2]) |>
    relocate(dataset) %>%
    split(f = .$dataset) |>
    map(
        .f = \(df) {
            df |>
                dplyr::filter(!sample %in% outliers_16s$outlier) |>
                select(-dataset) %>%
                select(., !caret::nearZeroVar(., names = TRUE)) # remove zero-var features
            }
    )

imap(
    t4f2_path_ls,
    \(x, i) write_csv(x, file.path("data", "processed", "16s", paste0("6_", i, "_pathways.csv")))
)

## name the 16S predicted pathways and compare with LC-MS putative pathways
### grab the convenient list of KEGG path names from the limma package
kegg_IDs <- limma::getKEGGPathwayNames() |> 
    rename(kegg_map = PathwayID, kegg_name = Description) |>
    mutate(kegg_ko = gsub("map", "ko", kegg_map)) |>
    relocate(kegg_name)

### make large table of LC-MS KEGG pathways, 16S KEGG pathways, and the KEGG pathID reference
seq_path_ls <- map(t4f2_path_ls, \(x) colnames(x)[-1]) |>
    unlist() |>
    unname() |>
    unique() %>%
    data.frame(kegg_ko = .) |>
    mutate(kegg_code = gsub("ko", "", kegg_ko))

lcms_path_ls <- lcms_pathways[names(t4f2_path_ls)] |>
    bind_rows() |>
    select(-sample) |>
    colnames() |>
    unique() |>
    data.frame() |>
    rename(kegg_name = 1) |>
    left_join(y = select(kegg_IDs, kegg_name, kegg_map), by = "kegg_name") |>
    mutate(kegg_code = gsub("map", "", kegg_map))

total_unique_paths <- full_join(
    x = lcms_path_ls,
    y = seq_path_ls,
    by = "kegg_code"
)

write_tsv(total_unique_paths, file.path("output", "tables", "6_KEGG_path_total.tsv"))

overlap_paths <- inner_join(
    x = lcms_path_ls,
    y = seq_path_ls,
    by = "kegg_code"
) |> 
    relocate(-kegg_code)

write_tsv(overlap_paths, file.path("output", "tables", "6_KEGG_path_overlaps.tsv"))

## perform sparse CCA by data set
set.seed(123)

cca_by_set <- function(input_x, i, input_z = t4f2_path_ls) {
    
    ### prepare data format
    df_x <- input_x |>
        semi_join(y = input_z[[i]], by = "sample") |>
        arrange(sample) |>
        column_to_rownames("sample")

    df_z <- input_z[[i]] |>
        semi_join(y = input_x, by = "sample") |>
        arrange(sample) |>
        column_to_rownames("sample")

    ### prepare the model
    permutation <- PMA::CCA.permute(
        x = df_x,
        z = df_z,
        typex = "standard",
        typez = "standard",
        nperms = 25,
        niter = 5,
        standardize = TRUE
    )

    ### save sparcification parameters
    penXtemp <- permutation$bestpenaltyx
    penZtemp <- permutation$bestpenaltyz

    ### conduct sparse CCA
    cca_res <- CCA(
        x = df_x,
        z = df_z,
        typex = "standard",
        typez = "standard",
        penaltyx = penXtemp,
        penaltyz = penZtemp,
        K = 3,
        niter = 15,
        standardize = TRUE
    )

    ### export canonical variates
    cca_res$u |> # LC-MS pathway data
        as.data.frame() |>
        rename(U1 = 1, U2 = 2, U3 = 3) |>
        mutate(variable = cca_res$xnames) |>
        relocate(variable) |>
        write_tsv(file.path("output", "tables", paste0("6_cca_", i, "_loadings_U_lcms.tsv")))
    
    cca_res$v |> # 16S pathway data
        as.data.frame() |>
        rename(V1 = 1, V2 = 2, V3 = 3) |>
        mutate(variable = cca_res$znames) |>
        relocate(variable) |>
        write_tsv(file.path("output", "tables", paste0("6_cca_", i, "_loadings_V_16s.tsv")))
    
    ### export correlations
    cca_res$cors %>%
        data.frame(correlations = .) |>
        write_csv(file.path("output", "tables", paste0("6_cca_", i, "_correlations.csv")))

    ### recast data along CV dimensions
    X_used <- as.matrix(df_x[, cca_res$xnames, drop = FALSE])
    Z_used <- as.matrix(df_z[, cca_res$znames, drop = FALSE])

    cca_scores <- as.data.frame(
        cbind(
            scale(X_used) %*% cca_res$u,
            scale(Z_used) %*% cca_res$v
        )
    ) |>
        rename(U1 = 1, U2 = 2, U3 = 3, V1 = 4, V2 = 5, V3 = 6) |>
        rownames_to_column("sample") |>
        as_tibble()

    write_csv(
        cca_scores,
        file.path("output", "tables", paste0("6_cca_", i, "_scores.csv"))
    )

    return(cca_scores)
}

cca_results <- imap(
    .x = lcms_pathways[names(t4f2_path_ls)],
    .f = cca_by_set
) |>
    bind_rows(.id = "dataset") |>
    tidyr::pivot_longer(cols = where(is.numeric), names_to = "cv", values_to = "scores") |>
    separate_wider_position(cv, c(cv_type = 1, cv_num = 1)) |>
    group_by(dataset, cv_num)

cca_key <- group_keys(cca_results)

cca_list <- group_split(cca_results)

## plot and export CCA results
### retrieve correlations
correlations <- list.files(path = file.path("output", "tables"),
    pattern = "^6_cca_(.*)\\_correlations.csv$",
    full.names = TRUE
) |>
    set_names(~ basename(.x) %>% sub("^6_cca_(.*)\\_correlations.csv$", "\\1", .)) |> 
    map(.f = \(x) read_csv(x) |> rownames_to_column("cv_num")) |>
    bind_rows(.id = "dataset") |>
    group_by(dataset, cv_num)

corr_key <- group_keys(correlations)

corr_list <- group_split(correlations)

### function for plotting CCA results
plot_cca <- function(df, i) {
    
    p <- df |>
        tidyr::pivot_wider(names_from = "cv_type", values_from = "scores") |>
        left_join(y = metadata, by = "sample") |>
        ggplot(aes(x = V, y = U)) +
            geom_point(aes(shape = culture_day,  color = exposure,), size = 1.75, alpha = 0.8, stroke = 1) +
            scale_color_manual(
                values = exposure_color_palette,
                drop = TRUE
            ) +
            scale_shape_manual(
                values = culture_day_shapes,
                drop = TRUE
            ) +
            theme_classic(12) +
            theme(
                legend.position = "none",
                plot.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "cm"),
                plot.title = element_text(size = 12)
            ) +
            coord_fixed(ratio = 1) +
            geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
            geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
            labs(
                x = paste0("16S variate ", unique(df$cv_num)),
                y = paste0("LC-MS variate ", unique(df$cv_num)),
                title = corr_list[[i]]$dataset,
                subtitle = paste0("corr. = ", round(corr_list[[i]]$correlations, 4))
            )
    
    return(p)
}

p_cca_ls <- imap(cca_list, plot_cca)

### prepare shared legend
legend_df <- expand.grid(
    exposure = levels(metadata$exposure),
    culture_day = levels(metadata$culture_day)
) |>
    dplyr::filter(culture_day != "d1") |>
    droplevels()

legend_plot <- ggplot(legend_df, aes(x = 1, y = 1)) +
    geom_point(aes(shape = culture_day, color = exposure), size = 2) +
    scale_color_manual(
        values = exposure_color_palette,
        drop = FALSE
    ) +
    scale_shape_manual(
        values = culture_day_shapes,
        drop = TRUE
    ) +
    theme_void() +
    theme(
        legend.position = "right",
        legend.title = element_text(face = "bold")
    )

shared_legend <- get_legend(legend_plot)

### assemble plot panel
pl <- cowplot::plot_grid(
    plotlist = p_cca_ls,
    labels = "AUTO",
    nrow = 3
)

pl_final <- plot_grid(
    pl, shared_legend,
    ncol = 2,
    rel_widths = c(1, 0.25)
)

ggsave(
    filename = file.path("output", "figures", "6_cca_pathways.png"),
    plot = pl_final,
    bg = "white",
    dpi = 300
)
