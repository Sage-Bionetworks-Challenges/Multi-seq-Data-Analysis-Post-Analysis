suppressPackageStartupMessages({
  library(parallel)
  library(data.table)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(tidyr)
  # library(writexl)
})

parser <- argparse::ArgumentParser()
parser$add_argument("-o", "--output_dir", help = "Output folder")
args <- parser$parse_args()

input_dir <- args$output_dir
output_dir <- args$output_dir

cluster_files <- list.files(input_dir, pattern = "*_comparisons.rds")
lfc_files <- list.files(input_dir, pattern = "*_lfc_results.rds")
umap_files <- list.files(input_dir, pattern = "*_umap.rds")

cluster_res <- mclapply(cluster_files, function(data) {
  readRDS(file.path(input_dir, data))
})

lfc_res <- mclapply(lfc_files, function(data) {
  readRDS(file.path(input_dir, data))
})

umap_res <- mclapply(umap_files, function(data) {
  readRDS(file.path(input_dir, data))
})

regex <- "(USF_biostat|zoradeng|moqri|BBKCS|DLS5|GOAL_LAB|LDExplore|Metformin-121|Baseline_MAGIC|Baseline_DeepImpute)_(ds\\d+[a-d]?)_(p\\d+k|p\\d+)"

table1 <- cluster_res %>% lapply(function(single_res) {
  map_df(single_res, ~map_df(.x, ~.x$cluster_metrics)) %>%
    select(query_name, reference_name, resolution, ari, nmi, n_clusters.query, n_clusters.ref) %>%
    tidyr::extract(query_name, into = c("team", "dataset", "proportion"), regex = regex, remove = FALSE) %>%
    select(team, dataset, proportion, everything())
}) %>% rbindlist()

table2 <- cluster_res %>% lapply(function(single_res) {
  map_df(single_res, ~map_df(.x, ~.x$knn_tibble)) %>%
  select(query_name, reference_name, barcode, shared_nn, union_nn, jaccard_shared) %>%
  group_by(query_name) %>%
  mutate(knn_jaccard_mean = mean(jaccard_shared)) %>%
  ungroup() %>%
  tidyr::extract(query_name, into = c("team", "dataset", "proportion"), regex = regex, remove = FALSE) %>%
  select(team, dataset, proportion, everything())
}) %>% rbindlist()

table3 <- cluster_res %>% lapply(function(single_res) {
  map_df(single_res, ~map_df(.x, ~.x$drop_outs)) %>%
  select(query_name, reference_name, tp, fp, tn, fn, f1) %>%
  tidyr::extract(query_name, into = c("team", "dataset", "proportion"), regex = regex, remove = FALSE) %>%
  select(team, dataset, proportion, everything())
}) %>% rbindlist()

table4 <- lfc_res %>% lapply(function(single_res) {
  combined_df <- bind_rows(
    df_10k <- single_res[which(names(single_res) =="lfc_10k")] %>% 
      rbindlist() %>%
      select(team, gene, avg_log2FC_1, p_val_1,  p_val_adj_1, avg_log2FC_10k, p_val_10k, p_val_adj_10k, comparsion_group) %>%
      pivot_longer(
        cols = -c(team, gene, comparsion_group), 
        names_to = c(".value", "resolution"), 
        names_pattern = "(.*)_(.*)"),
    df_20k <- single_res[which(names(single_res) =="lfc_20k")] %>% 
      rbindlist() %>%
      select(team, gene, avg_log2FC_1, p_val_1,  p_val_adj_1, avg_log2FC_20k, p_val_20k, p_val_adj_20k, comparsion_group) %>%
      pivot_longer(
        cols = -c(team, gene, comparsion_group), 
        names_to = c(".value", "resolution"), 
        names_pattern = "(.*)_(.*)"),
    df_50k <- single_res[which(names(single_res) =="lfc_50k")] %>% 
      rbindlist() %>%
      select(team, gene, avg_log2FC_1, p_val_1,  p_val_adj_1, avg_log2FC_50k, p_val_50k, p_val_adj_50k, comparsion_group) %>%
      pivot_longer(
        cols = -c(team, gene, comparsion_group), 
        names_to = c(".value", "resolution"), 
        names_pattern = "(.*)_(.*)")
  )
}) %>% rbindlist()

## UMAP
ds_reads_props <- c("p10k", "p20k", "p50k")
query_umap_df <- umap_res %>% lapply(function(single_res) {
  map_df(single_res, ~map_df(.x, ~bind_rows(.x))) %>%
    select(query, cell_barcode, UMAP_1, UMAP_2, nCount_RNA) %>% 
    rename(query_name = query) %>%
    tidyr::extract(query_name, into = c("team", "dataset", "proportion"), regex = regex, remove = FALSE) %>%
    select(team, dataset, proportion, everything()) %>%
    mutate(downsampled_by = case_when(proportion %in% ds_reads_props ~ "reads",
                                     TRUE ~ "cells"))
}) %>% rbindlist()

# get ref's umap info
# datasets <- c("ds1a", "ds1b", "ds1c", "ds1d", "ds2", "ds3")
# ref_umap_df <- pbmclapply(datasets, function(dataset) {
  
#   f <- str_glue("{dataset}_p1")
#   reference <- readRDS(str_glue("seurat_objects/{f}.rds"))

#   # get library size factor
#   sf_df <- FetchData(reference, vars = "nCount_RNA")
#   # create umap info table
#   reference[["umap"]]@cell.embeddings %>% 
#     as.data.frame() %>%
#     tibble::rownames_to_column("cell_barcode") %>% 
#     mutate(team = "no_imp",
#            dataset = dataset,
#            proportion = "p1",
#            query_name = reference@project.name,
#            nCount_RNA = sf_df$nCount_RNA[match(cell_barcode, rownames(sf_df))],
#            downsampled_by = "ref") %>%
#     select(colnames(query_umap_df))
#   }, mc.cores = length(datasets)) %>% rbindlist()

table5 <- query_umap_df
sheets <- list(
  "cluster_metrics" = table1,
  "knn_jaccard" = table2,
  "drop_outs" = table3,
  "lfc" = table4,
  "umap" = table5
)

saveRDS(sheets, file.path(output_dir, "downstream_results_summary.rds"))
# write_xlsx(sheets, file.path(output_dir, "downstream_results_summary.xlsx"))

system(str_glue("synapse store {file.path(output_dir, 'downstream_results_summary.rds')} --id syn52965907"))
# system(str_glue("synapse store {file.path(output_dir, 'downstream_results_summary.xlsx'} --id syn52965908"))
