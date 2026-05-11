# 13) Identifying age-dependent, exposure-differentiating pathways
## organize metadata for joining disparate results
meta_sets <- mutate(
    metadata,
    dataset = paste(cohort, dose, age_week, sep = "_")
) |>
    select(dataset, cohort, lcms_run) |>
    distinct()

annotation_key <- read_tsv(file.path("data", "processed", "lcms", "5_putative_annotation_key.tsv"))

## grab pathways that vary by age
lcms_de_inocula <- read_csv(file.path("output", "tables", "12_lcms_DE_paths_by_age_inocula.csv")) |>
    inner_join(
        y = select(meta_sets, lcms_run, cohort) |> distinct(),
        by = "lcms_run",
        relationship = "many-to-many"
    )

## grab pathways that vary by exposure
sig_path_annotations <- read_tsv(file.path("output", "tables", "5_lcms_sig_paths_annotated.tsv"))

## intersect age- and exposure-differential pathways
age_exp_paths <- sig_path_annotations |>
    select(dataset, culture_day, pathwayID = pathway, feature_counts) |>
    left_join(
        y = meta_sets,
        by = "dataset"
    ) |>
    inner_join( 
        y = lcms_de_inocula,
        by = c("cohort", "pathwayID" = "pathway"),
        relationship = "many-to-many"
    ) |>
    left_join(
        y = select(annotation_key, pathwayID, pathway),
        by = c("pathwayID"),
        relationship = "many-to-many"
    ) |>
    relocate(rep_features = feature_counts, pathway, .after = last_col()) |>
    arrange(dataset, pathwayID, culture_day) |>
    drop_na() |>
    select(-cohort, -p_age, -starts_with("lcms_run")) |>
    distinct()

## export summaries of the intersections
write_tsv(
    age_exp_paths,
    file.path("output", "tables", "13_age_exp_pathways.tsv")
)

summary_by_set <- age_exp_paths |>
    group_by(dataset, culture_day, direction) |>
    summarize(path_count = n(), paths = paste(pathwayID, collapse = "; "), .groups = "drop_last") |>
    ungroup()

write_tsv(summary_by_set, file.path("output", "tables", "13_age_exp_summary_by_set.tsv"))

age_exp_paths |>
    group_by(pathwayID, pathway) |>
    summarize(set_count = n(), sets = paste(dataset, collapse = "; "), .groups = "drop_last") |>
    arrange(desc(set_count)) |>
    write_tsv(file.path("output", "tables", "13_age_exp_summary_by_pathway.tsv"))

## Plotting the pathway counts and direction by $dataset and $culture_day
### organize data as list
summary_by_set_meta <- summary_by_set |>
    separate_wider_delim(cols = dataset, delim = "_", names = c("cohort", "dose", "age_week")) |>
    mutate(
        culture_day = factor(culture_day, levels = levels(metadata$culture_day)),
        culture_day = fct_drop(culture_day),
        dose = factor(dose, levels = c("low", "high")),
        age_week = factor(age_week, levels = levels(metadata$age_week)),
        cohort = factor(cohort)
    ) |>
    group_by(age_week, cohort, .drop = FALSE)

summary_by_set_key <- group_keys(summary_by_set_meta)

list_by_set <- group_split(summary_by_set_meta)

### prepare list of plots
p_summary_by_set <- imap(
    .x = list_by_set,
    .f = \(set, i) {

        title_i <- paste(
            summary_by_set_key$cohort[i],
            summary_by_set_key$age_week[i]
        )

        if (nrow(set) == 0) {
            ggplot() +
                theme_void() +
                labs(subtitle = title_i) +
                annotate("text", x = 0.5, y = 0.5,
                    label = "No age- and exposure-dependent\npathways",
                    hjust = 0.5, vjust = 0.5, size = 3) +
                theme(plot.subtitle = element_text(size = 8.5, margin = margin(l = 40, t = 6)))

        } else {

        ggplot(set, aes(x = culture_day, y = dose)) +
            geom_point(
                aes(size = path_count, fill = direction),
                shape = 21,
                position = position_dodge(width = 0.4)
            ) +
            scale_x_discrete(
                name = "culture day",
                drop = FALSE
            ) +
            scale_y_discrete(
                drop = FALSE
            ) +
            ggplot2::scale_size_continuous(
                name = "Number of pathways",
                limits = c(1, 45)
            ) +
            scale_fill_manual(
                name = "Trend across inocula",
                values = c(
                    increasing_with_age = "#1b9e77",
                    decreasing_with_age = "#d95f02"
                )
            ) +
            theme_cowplot(10) +
            background_grid() +
            theme(axis.text.x  = element_text(angle = 0, hjust = 0.5), legend.position = "none") +
            labs(subtitle = title_i)
        }
    }
)

### prepare legend
legend_df <- data.frame(
    path_count = c(5, 40),
    direction = c("increasing_with_age", "decreasing_with_age"),
    culture_day = factor(c("d1", "d2"), levels = levels(summary_by_set_meta$culture_day)),
    dose = factor(c("low", "high"), levels = c("low", "high"))
)

legend_plot <- ggplot(legend_df, aes(x = culture_day, y = dose)) +
    geom_point(
        aes(size = path_count, fill = direction),
        shape = 21,
        position = position_dodge(width = 0.4)
    ) +
    scale_size_continuous(name = "Pathway\ncount") +
    scale_fill_manual(
        name = "Trend across\ninocula",
        values = c(
            increasing_with_age = "#1b9e77",
            decreasing_with_age = "#d95f02"
        )
    ) +
    theme_cowplot() +
    theme(
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 8)
    )

lgn <- get_legend(legend_plot)

### arrange the panel and export
p_grid <- plot_grid(
    plotlist = p_summary_by_set,
    nrow = 5,
    labels = "AUTO"
)

panel <- plot_grid(
    p_grid,
    ggdraw(lgn),
    align = "h",
    nrow = 1, ncol = 2,
    rel_widths = c(1, 0.3)
)

ggsave(
    filename = file.path("output", "figures", "13_age_and_exposure.png"),
    plot = panel,
    dpi = 300,
    bg = "white",
    height = 7,
    width = 6.5
)
