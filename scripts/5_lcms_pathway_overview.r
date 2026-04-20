# 5) LC-MS putative pathways analyzed by limma
if(!exists("lcms_pathways", mode = "list")) {
    lcms_pathways <- list.files(path = file.path("data", "processed", "lcms"),
        pattern = "^4_",
        full.names = TRUE
    ) |>
    set_names(~ basename(.x) %>% sub("^1_(.*)\\.csv$", "\\1", .)) |> 
    map(.f = read_csv)
}

## limma workflow
