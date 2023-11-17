suppressPackageStartupMessages({
  library(parallel)
  library(data.table)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(tidyr)
  library(writexl)
})

input_dir <- "/input"
output_dir <- "/output"

cluster_files <- list.files(input_dir, pattern = "*_comparsion.rds")
lfc_files <- list.files(input_dir, pattern = "*_lfc_results.rds")

cluster_res <- mclapply(cluster_files, function(data) {
  readRDS(file.path(input_dir, data))
})

lfc_res <- mclapply(lfc_files, function(data) {
  readRDS(file.path(data))
})

regex <- "(USF_biostat|zoradeng|moqri|BBKCS|DLS5|GOAL_LAB|LDExplore|Metformin-121)_(ds\\d+[a-d]?)_(p\\d+k|p\\d+)"

table1 <- cluster_res %>% lapply(function(single_res) {
  map_df(single_res, ~map_df(.x, ~.x$cluster_metrics)) %>%
    select(query_name, reference_name, resolution, ari, nmi, n_clusters.query, n_clusters.ref) %>%
    extract(query_name, into = c("team", "dataset", "proportion"), regex = regex, remove = FALSE) %>%
    select(team, dataset, proportion, everything())
}) %>% rbindlist()

table2 <- cluster_res %>% lapply(function(single_res) {
  map_df(single_res, ~map_df(.x, ~.x$knn_tibble)) %>%
  select(query_name, reference_name, barcode, shared_nn, union_nn, jaccard_shared) %>%
  group_by(query_name) %>%
  mutate(knn_jaccard_mean = mean(jaccard_shared)) %>%
  ungroup() %>%
  extract(query_name, into = c("team", "dataset", "proportion"), regex = regex, remove = FALSE) %>%
  select(team, dataset, proportion, everything())
}) %>% rbindlist()

table3 <- cluster_res %>% lapply(function(single_res) {
  map_df(single_res, ~map_df(.x, ~.x$drop_outs)) %>%
  select(query_name, reference_name, tp, fp, tn, fn, f1) %>%
  extract(query_name, into = c("team", "dataset", "proportion"), regex = regex, remove = FALSE) %>%
  select(team, dataset, proportion, everything())
}) %>% rbindlist()

table4 <- lfc_res %>% lapply(function(single_res) {
  combined_df <- bind_rows(
    df_10k <- single_res[which(names(single_res) =="lfc_10k")] %>% 
      rbindlist() %>%
      select(team, avg_log2FC_1, p_val_1,  p_val_adj_1, avg_log2FC_10k, p_val_10k, p_val_adj_10k, comparsion_group) %>%
      pivot_longer(
        cols = -c(team, comparsion_group), 
        names_to = c(".value", "resolution"), 
        names_pattern = "(.*)_(.*)"),
    df_20k <- single_res[which(names(single_res) =="lfc_20k")] %>% 
      rbindlist() %>%
      select(team, avg_log2FC_1, p_val_1,  p_val_adj_1, avg_log2FC_20k, p_val_20k, p_val_adj_20k, comparsion_group) %>%
      pivot_longer(
        cols = -c(team, comparsion_group), 
        names_to = c(".value", "resolution"), 
        names_pattern = "(.*)_(.*)"),
    df_50k <- single_res[which(names(single_res) =="lfc_50k")] %>% 
      rbindlist() %>%
      select(team, avg_log2FC_1, p_val_1,  p_val_adj_1, avg_log2FC_50k, p_val_50k, p_val_adj_50k, comparsion_group) %>%
      pivot_longer(
        cols = -c(team, comparsion_group), 
        names_to = c(".value", "resolution"), 
        names_pattern = "(.*)_(.*)")
  )
}) %>% rbindlist()

sheets <- list(
  "cluster_metrics" = table1,
  "knn_jaccard" = table2,
  "drop_outs" = table3,
  "lfc" = table4
)

saveRDS(sheets, file.path(output_dir, "downstream_results_summary.rds"))
write_xlsx(sheets[-4], file.path(output_dir, "downstream_results_summary.xlsx"))

system(str_glue("synapse store {file.path(output_dir, 'downstream_results_summary.rds'} --id syn52965907"))
system(str_glue("synapse store {file.path(output_dir, 'downstream_results_summary.xlsx'} --id syn52965908"))
