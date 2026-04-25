# 6) CCA on LC-MS log2 peak intensities and 16S predicted pathway counts
if(!exists("lcms_pathways", mode = "list")) {
    lcms_pathways <- list.files(path = file.path("data", "processed", "lcms"),
        pattern = "^4_",
        full.names = TRUE
    ) |>
    set_names(~ basename(.x) %>% sub("^4_(.*)\\.tsv$", "\\1", .)) |> 
    map(.f = read_tsv)
}

t4f2_pathways <- read_csv(file.path("data", "processed", "16s", "sequencing_Tax4Fun2_pathways.csv"))

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

## name the 16S predicted pathways and compare with LC-MS putative pathways
### grab the convenient list of KEGG path names from the limma package
kegg_names <- limma::getKEGGPathwayNames() |> 
    rename(kegg_pathID = PathwayID, kegg_name = Description)

### make large table of LC-MS KEGG pathways, 16S KEGG pathways, and the KEGG pathID reference

## perform sparse CCA by data set
set.seed(123)

cca_by_set <- function(input_x, i, input_z = t4f2_path_ls) {

    ### prepare data format
    df_x <- input_x |>
        semi_join(y = input_z[[i]], by = "sample") |>
        column_to_rownames("sample")

    df_z <- input_z[[i]] |>
        semi_join(y = input_x, by = "sample") |>
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
    cca_scores <- as.data.frame(
        cbind(
            as.matrix(select(df_x, cca_res$xnames) %>% scale()) %*% cca_res$u,
            as.matrix(select(df_z, cca_res$znames) %>% scale()) %*% cca_res$v
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

cca_list <- imap(
    .x = lcms_pathways[names(t4f2_path_ls)],
    .f = cca_by_set
)

## plot CCA
