# 12) Comparing putative metabolic pathways across inocula from different ages
if(!exists("lcms_pathways", mode = "list")) {
    lcms_pathways <- list.files(path = file.path("data", "processed", "lcms"),
        pattern = "^4_(BPS|DINP)_",
        full.names = TRUE
    ) |>
    set_names(~ basename(.x) %>% sub("^4_(.*)\\.tsv$", "\\1", .)) |> 
    map(.f = read_tsv)
}

## subset data to BPS and DINP high-exposure LC-MS runs
meta_high_ex_d0 <- metadata |>
    dplyr::filter(
        grepl("high-bps|high-dinp", lcms_run),
        exposure == "inoculum"
    ) |>
    mutate(dataset = paste(cohort, dose, age_week, sep = "_")) |>
    droplevels()

paths_high_ex_d0 <- lcms_pathways[names(lcms_pathways) %in% meta_high_ex_d0$dataset] |>
    bind_rows() |>
    inner_join(y = select(meta_high_ex_d0, sample, lcms_run), by = "sample") |>
    group_by(lcms_run)

lcms_run_key <- group_keys(paths_high_ex_d0)

paths_list_d0 <- paths_high_ex_d0 |> 
    group_split(.keep = FALSE) |>
    map(.f = \(x) {x %>% select(where(~ any(!is.na(.))))})

## set up linear mixed-effects model per pathway
fit_pathway <- function(y, md) {

    md$path <- as.numeric(y)

    model <- lmerTest::lmer(
        path ~ age_week + (1 | replicate), 
        data = md,
        REML = FALSE
    )

    ### retrieve p-value for the pathway
    p_age <- as.numeric(anova(model)["Pr(>F)"])

    ### estimated marginal means by age
    em <- emmeans(model, ~ age_week, mode = "satterthwaite")
    em_df <- as.data.frame(em)[, c("age_week", "emmean")]

    list(
        p_age  = p_age,
        emmean = em_df
    )
}

## set up function for summarizing directionality from emmeans
age_numeric <- c(
    "wk07" = 7,
    "wk12" = 12,
    "wk17" = 17,
    "wk22" = 22
)

summarize_direction <- function(em_df) {
    x <- age_numeric[as.character(em_df$age_week)]
    y <- em_df$emmean

    slope <- coef(lm(y ~ x))[2]

    direction <- dplyr::case_when(
        slope > 0 ~ "increasing_with_age",
        slope < 0 ~ "decreasing_with_age",
        TRUE ~ "flat"
    )

    list(slope = slope, direction = direction)
}

## apply function to the datasets
lme_res_d0 <- map(
    .x = paths_list_d0,
    .f = \(x) {

        meta <- semi_join(
            meta_high_ex_d0,
            x,
            by = "sample")
        
        paths <- select(x, starts_with("map")) |>
            as.matrix()
        
        ### run fit per path across the data matrix
        path_fits <- apply(
            paths,
            2,
            fit_pathway,
            md = meta
        )

        ### compute FDR
        age_pval <- sapply(path_fits, '[[', "p_age")
        age_fdr <- p.adjust(age_pval, method = "fdr")

        ### extract emmeans info from list formats
        dir_list <- lapply(path_fits, \(res) summarize_direction(res$emmean))
        slopes <- sapply(dir_list, '[[', "slope")
        dirs <- sapply(dir_list, '[[', "direction") 

        ### tabulate all results for listing
        pathway_results <- tibble(
            pathway = colnames(paths),
            p_age = age_pval,
            FDR_age = age_fdr,
            slope = slopes,
            direction = dirs
        )
    }
)

## export pathways which differ significantly between the inocula of mice of different ages
names(lme_res_d0)  <- lcms_run_key$lcms_run
lme_res_d0 |>
    bind_rows(.id = "lcms_run") |>
    dplyr::filter(FDR_age < .05) |>
    write_csv(file.path("output", "tables", "12_lcms_DE_paths_by_age_inocula.csv"))
