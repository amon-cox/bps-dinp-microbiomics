# 7) Taxonomic profiling
## grab 16S sequencing data
rel_ab_processed <- list.files(path = file.path("data", "processed", "16s"),
    pattern = "^1_",
    full.names = TRUE
) |>
    set_names(~ basename(.x) %>% sub("^1_(.*)\\.csv$", "\\1", .)) |> 
    map(.f = read_csv)

## find average abundance of ASVs per group per dataset
avg_rel_ab <- map(
    .x = rel_ab_processed,
    .f = \(df) {
        left_join(
            x = df,
            y = select(metadata, sample, exposure, culture_day),
            by = "sample"
        ) |>
        pivot_longer(
            cols = where(is.numeric),
            names_to = "ASV",
            values_to = "rel_abnd"
        ) |>
        summarize(
            asv_avg_abnd = mean(rel_abnd),
            .by = c(ASV, exposure, culture_day)
        ) |>
        arrange(culture_day, exposure, desc(asv_avg_abnd))
    }
) |>
    bind_rows(.id = "dataset")

## get taxonomic labels for included ASVs
taxons <- read_csv(file.path("data", "original", "16s", "sequencing_feature_metadata.csv"))

taxon_labels <- taxons |>
    dplyr::filter(ASV %in% unique(avg_rel_ab$ASV)) %>%
    mutate(
        across(Domain:Species, ~ ifelse(. == "__", NA, .)),
        across(Domain:Species, ~ ifelse(grepl("uncultured", .), NA, .)),
        best_taxon = coalesce(Genus, Family, Order, Class, Phylum, Domain)
    ) |>
    select(ASV, best_taxon)

## filter to top 20 most abundant ASVs in the inoculum, across datasets
top_asvs <- avg_rel_ab |>
    dplyr::filter(exposure == "inoculum") |>
    group_by(dataset) |>
    arrange(desc(asv_avg_abnd)) |>
    slice_head(n = 20) |>
    ungroup() |>
    select(ASV)

asvs_to_plot <- avg_rel_ab |>
    left_join(y = taxon_labels, by = "ASV") |>
    dplyr::filter(ASV %in% top_asvs$ASV) |>
    droplevels() |>
    mutate(
        ASV_taxon = paste0("[", ASV, "] ", best_taxon),
        asv_avg_abnd_per = round(asv_avg_abnd * 100, 2),
        ASV_taxon = forcats::fct_reorder(ASV_taxon, asv_avg_abnd)
    )

## plot average relative abundances as heatmap
p_asvs_avg <- ggplot(
    asvs_to_plot,
    aes(x = culture_day, y = ASV_taxon, fill = asv_avg_abnd_per, group = dataset)
) +
    geom_tile() +
    geom_text(aes(label = asv_avg_abnd_per), color = "white", size = 2) +
    facet_wrap(vars(dataset, exposure), scales = "free_x") +
    scale_fill_continuous(palette = 'viridis', limits = c(0, 50)) +
    theme_cowplot() +
    theme(
        strip.text = element_text(face = "bold", hjust = 0.5, size = 8),
        strip.background = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        panel.spacing.x = unit(0.15, "cm"),
        plot.title = element_text(hjust = 1),
        axis.text.y = element_text(hjust = 0, size = 6),
        axis.text.x = element_text(size = 6),
        axis.line = element_blank(),
        axis.ticks = element_blank()
    ) +
    labs(
        x = element_blank(),
        y = element_blank(),
        fill = paste0("Average (%)\nRelative\nAbundance")
    )

ggsave(
    filename = file.path("output", "figures", "7_average_rel_ab.png"),
    plot = p_asvs_avg,
    bg = "white",
    dpi = 300,
    height = 8,
    width = 6
)
