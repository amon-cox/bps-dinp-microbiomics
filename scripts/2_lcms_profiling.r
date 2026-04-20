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

imap(.x = lcms_distances, .f = permanova_by_set)

## plot each cohort and dose by week, then arrange into panels
### prepare plotting function
pcoa_by_set <- function(ds, ds_name) {

    parts <- strsplit(ds_name, "_") |> unlist()
    cohort <- parts[1]
    dose <- parts[2]
    age_week <- parts[3]

    #### multidimensional scaling, i.e. principal coordinates analysis
    pcoa <- cmdscale(d = ds, k = 2) |>
        as.data.frame() |>
        rownames_to_column("sample") |>
        left_join(y = metadata, by = "sample")

    #### plotting
    p <- ggplot(pcoa, aes(x = V1, y = V2, fill = exposure)) +
        geom_point(aes(shape = culture_day,  color = exposure,), size = 1.75, alpha = 0.8, stroke = 1) +
        scale_color_manual(
            aesthetics = c("color", "fill"),
            values = exposure_color_palette,
            drop = TRUE
        ) +
        scale_shape_manual(
            values = culture_day_shapes,
            drop = TRUE
        ) +
        stat_ellipse(geom = "polygon", alpha = 0.2) +
        theme_classic(12) +
        theme(
            legend.position = "none",
            plot.margin = margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm")
        ) + 
        coord_fixed(
            ratio = 1,
            xlim = range(pcoa[ , c("V1", "V2")]) * 1.5,
            ylim = range(pcoa[ , c("V1", "V2")]) * 1.5
        ) +
        labs(x = "component 1", y = "component 2", title = ds_name)
    
    #### save info for setting up panel later
    tbl <- tibble(
        dataset = ds_name,
        cohort = cohort,
        dose = dose,
        age_week = age_week,
        plot = list(p)
    )
}

### store all plotting info
pcoa_tbl <- map2_dfr(
    .x = lcms_distances,
    .y = names(lcms_distances),
    .f = pcoa_by_set
)

### arrange panels by cohort x dose
#### 100uM BPS panel
panel_bps <- pcoa_tbl |>
    filter(cohort == "BPS", dose == "high") |>
    arrange(age_week) |>
    pull(plot) |>
    (\(.x) {
        legend <- get_legend(.x[[1]] + theme(
            legend.position = "right",
            legend.title = element_text(face = "bold")
        ))

        pl <- plot_grid(
            plotlist = .x,
            labels = c("A", "B", "C", "D", "E", ""),
            nrow = 2
        )
        
        pl_final <- pl + draw_grob(
            legend, 
            x = 4.6/6, y = 1/6, width = 0.2, height = 0.2
        )
    })()

ggsave(
    filename = file.path("output", "figures", "2_pcoa_high-BPS.png"),
    plot = panel_bps,
    bg = "white",
    height = 5.75,
    width = 8,
    dpi = 300
)

#### 100uM DINP panel
panel_dinp <- pcoa_tbl |>
    filter(cohort == "DINP", dose == "high") |>
    arrange(age_week) |>
    pull(plot) |>
    (\(.x) {
        legend <- get_legend(.x[[1]] + theme(
            legend.position = "right",
            legend.title = element_text(face = "bold")
        ))

        pl <- plot_grid(
            plotlist = .x,
            labels = c("A", "B", "C", "D", "E", ""),
            nrow = 2
        )
        
        pl_final <- pl + draw_grob(
            legend, 
            x = 4.6/6, y = 1/6, width = 0.2, height = 0.2
        )
    })()

ggsave(
    filename = file.path("output", "figures", "2_pcoa_high-DINP.png"),
    plot = panel_dinp,
    bg = "white",
    height = 5.75,
    width = 8,
    dpi = 300
)

#### 10uM panel
#### dummy plot for shared legend
legend_df <- expand.grid(
    exposure = levels(metadata$exposure),
    culture_day = levels(metadata$culture_day)
) |>
    dplyr::filter(culture_day != "d1") |>
    droplevels()

legend_plot <- ggplot(legend_df, aes(x = 1, y = 1, fill = exposure)) +
    geom_point(aes(shape = culture_day, color = exposure), size = 2) +
    scale_color_manual(
        aesthetics = c("colour", "fill"),
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

#### plot the 10uM panel
panel_low <- pcoa_tbl |>
    filter(dose == "low") |>
    arrange(cohort, age_week) |>
    pull(plot) |>
    (\(.x) {
        pl <- plot_grid(
            plotlist = .x,
            labels = c("A", "B", "C", "D", "E", ""),
            nrow = 2
        )
        
        pl_final <- pl + draw_grob(
            shared_legend, 
            x = 4.6/6, y = 1/6, width = 0.2, height = 0.2
        )
    })()

ggsave(
    filename = file.path("output", "figures", "2_pcoa_low-select.png"),
    plot = panel_low,
    bg = "white",
    height = 5.75,
    width = 8,
    dpi = 300
)
