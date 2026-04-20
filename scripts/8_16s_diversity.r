# 8) Diversity analyses for 16S sequencing data
if(!exists("rel_ab_processed", mode = "list")) {
    rel_ab_processed <- list.files(path = file.path("data", "processed", "16s"),
        pattern = "^1_",
        full.names = TRUE
    ) |>
        set_names(~ basename(.x) %>% sub("^1_(.*)\\.csv$", "\\1", .)) |> 
        map(.f = read_csv)
}

## alpha diversity
set.seed(123)

alpha_div <- imap(
    .x = rel_ab_processed,
    .f = \(df, nm) {
        ### calculate alpha-diversity scores
        alpha <- df |>
            rowwise() |>
            mutate(
                shannon = diversity(c_across(-sample), index = "shannon"),
                invsimpson = diversity(c_across(-sample), index = "invsimpson")
            ) |>
            ungroup() |>
            left_join(y = select(metadata, sample, exposure, culture_day), by = "sample") |>
            mutate(
                group = paste(exposure, culture_day, sep = "_"),
                group = factor(group)
            ) |>
            select(!starts_with("ASV"))
        
        write_csv(
            x = alpha,
            file = file.path("output", "tables", paste0("8_16s_alpha_scores_", nm, ".csv"))
        )

        ### Kruskal-Wallis test
        kw_test <- data.frame(
            shannon = kruskal.test(shannon ~ group, data = alpha)$p.value,
            invsimpson = kruskal.test(invsimpson ~ group, data = alpha)$p.value
        )

        ### post-hoc test with Dunn's test with Bonferroni correction
        shannon <- dunnTest(shannon ~ group, data = alpha, method = "bonferroni")$res
        invsimpson <- dunnTest(invsimpson ~ group, data = alpha, method = "bonferroni")$res

        dunn_test <- list(shannon, invsimpson) |>
            bind_rows(.id = "index")
        
        write_csv(
            x = dunn_test,
            file = file.path("output", "tables", paste0("8_16s_alphaDiv_DunnTest_", nm, ".csv"))
        )

        ### reorganize diversity info
        alpha_info <- pivot_longer(
                alpha,
                cols = c(shannon, invsimpson),
                names_to = "index",
                values_to = "score"
            ) |>
            mutate(
                group = forcats::fct_relevel(
                    group, 
                    "inoculum_d0", "control_d2", "control_d7"),
                index = forcats::fct_relevel(index, "shannon")
            )

        ### violin plot of alpha diversity scores, returned to the list
        p_alpha <- alpha_info |>
            ggplot(aes(x = group, y = score)) +
                geom_violin(
                    aes(fill = exposure, color = exposure),
                    alpha = 0.75, linewidth = 0.75, trim = FALSE
                ) +
                scale_color_manual(
                    aesthetics = c("color", "fill"),
                    values = exposure_color_palette,
                    drop = TRUE
                ) +
                theme_half_open(10) +
                theme(
                    legend.position = "none",
                    plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "in"),
                    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                    plot.title = element_text(size = 10)
                ) +
                facet_wrap(~ index, scales = "fixed") +
                scale_y_continuous(limits = c(0, NA), expand = c(0,0)) +
                labs(x = element_blank(), y = element_blank(), title = nm)
    }
)

## beta diversity
beta_div <- imap(
    .x = rel_ab_processed,
    .f = \(df, nm) {
        ### calculate dissimilarity matrix
        bray_dism <- df |>
            column_to_rownames("sample") |>
            vegdist(method = "bray") |>
            as.matrix()
        
        ### Principal Coordinates Analysis on dissimilarity matrix
        pcoa <- cmdscale(bray_dism, k = 2) |>
            as.data.frame() |>
            rownames_to_column("sample") |>
            left_join(select(metadata, sample, exposure, culture_day), by = "sample") |>
            mutate(
                group = paste(exposure, culture_day, sep = "_"),
                group = factor(group),
                group = forcats::fct_relevel(
                    group, 
                    "inoculum_d0", "control_d2", "control_d7")
            )
        
        ### conduct PERMANOVA by group for dissimilarities matrix
        permanova <- adonis2(bray_dism ~ group, data = pcoa)
        pairwise <- pairwiseAdonis::pairwise.adonis2(bray_dism ~ group, data = pcoa)

        write_csv(
            permanova,
            file.path("output", "tables", paste0("8_16s_betaDiv_PERMANOVA_", nm, ".csv"))
        )

        capture.output(
            pairwise,
            file = file.path("output", "tables", paste0("8_16s_betaDiv_pairwise_", nm, ".txt"))
        )

        ### prepare plots to return as list
        p_beta <- pcoa |>
            ggplot(aes(x = V1, y = V2, fill = exposure)) +
                geom_point(aes(color = exposure, shape = culture_day), size = 2, alpha = 0.65, stroke = 1.25) +
                scale_color_manual(
                    aesthetics = c("color", "fill"),
                    values = exposure_color_palette,
                    drop = TRUE
                ) +
                scale_shape_manual(values = c(16, 15, 3)) +
                theme_half_open(10) +
                theme(
                    legend.position = "none",
                    plot.title = element_text(size = 10)
                ) +
                coord_fixed() +
                stat_ellipse(geom = "polygon", level = 0.95, alpha = 0.15) +
                labs(x = "component 1", y = "component 2", title = nm, shape = "culture\nday")
    }
)

## arrange diversity plots into panel figure
### alpha-div column
pl_alpha <- plot_grid(
    plotlist = alpha_div,
    align = "hv",
    ncol = 1,
    labels = c("A", "B", "C")
)

### beta-div column
pl_beta <- plot_grid(
    plotlist = beta_div,
    align = "hv",
    ncol = 1,
    labels = c("D", "E", "F")
)

### shared legend
legend_df <- expand.grid(
    exposure = factor(c("inoculum", "control", "bps", "dinp"),
        levels = c("inoculum", "control", "bps", "dinp")),
    culture_day = factor(c("d0", "d2", "d7"),
        levels = c("d0", "d2", "d7"))
)

legend_plot <- ggplot(legend_df, aes(x = 1, y = 1, fill = exposure)) +
    geom_point(aes(shape = culture_day, color = exposure), size = 2) +
    scale_color_manual(
        aesthetics = c("color", "fill"),
        values = exposure_color_palette,
        drop = TRUE
    ) +
    scale_shape_manual(
        values = culture_day_shapes,
        drop = TRUE
    ) +
    stat_ellipse(geom = "polygon", level = 0.95, alpha = 0.15) +
    theme_void() +
    theme(
        legend.position = "right",
        legend.title = element_text(face = "bold", size = 10)
    )

legend <- get_legend(legend_plot)

### final figure
pl_div <- plot_grid(
    pl_alpha,
    pl_beta,
    legend,
    ncol = 3,
    rel_widths = c(1, 1, 0.25)
)

ggsave(
    filename = file.path("output", "figures", "8_diversity_scores.png"),
    plot = pl_div,
    bg = "white",
    dpi = 300
)
