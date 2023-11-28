# the script is run on the ec2-instance
# Rscript downstream_task1.R --ds_dir /challenge-data/downsampled/scrna/ \
#  --raw_dir /challenge \
#  -g /input/scrna_goldstandard.rds \
#  -f /input/DLS5/ \
#  --output_dir /evals/ \
#  --workers 20 \
#  --cores 6 
# args <- list()
# args$ds_dir = "/challenge-data/downsampled/scrna/"
# args$raw_dir = "/challenge"
# args$output_dir = "/evals"
# args$submission_dir = "/input/zoradeng/"
# args$goldstandard = "/input/scrna_goldstandard.rds"
# args$workers = 20
# args$cores = 6

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(data.table)
  library(magrittr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(aricode)
  library(parallel)
  library(pbmcapply)
  library(future)
  library(Seurat)
})

# Set up arguments
parser <- argparse::ArgumentParser()
parser$add_argument("--ds_dir", help = "Downsampled data folder")
parser$add_argument("--raw_dir", help = "Raw data folder")
parser$add_argument("-f", "--submission_dir", help = "Submission output folder")
parser$add_argument("-g", "--goldstandard", help = "Goldstandard file")
parser$add_argument("-o", "--output_dir", help = "Output folder")
parser$add_argument("--workers", help = "Workers for parallelization")
parser$add_argument("--cores", help = "Cores for parallelization")

args <- parser$parse_args()

plan("multicore", workers = as.numeric(args$workers)) # set up seurat funcs parallelization
cores <- as.numeric(args$cores)

source("downstream_funcs.R")
options(readr.show_col_types = FALSE)
options(future.seed = TRUE)
options(future.globals.maxSize = 3000 * 1024^2)  # increase to 3G
options(ignore.interactive = TRUE)

## Set up data paths
raw_dir <- args$raw_dir
ds_dir <-  args$ds_dir
data_dir <- args$submission_dir # submission_dir is named as team name
data_dir = args$ds_dir
team <- basename(data_dir) # get team name
team <- "no_imp"
output_dir <- args$output_dir
print(team)
# Set up datasets info
datasets <- c("ds1a", "ds1b", "ds1c", "ds1d", "ds2", "ds3")
ds_reads_props <- c("p10k", "p20k", "p50k")
ds_cells_props <- c("p00625", "p0125", "p025")

# Read valid cells
valid_cell_path <- file.path(output_dir, "valid_cells.rds")
if (!file.exists(valid_cell_path)) {
  valid_cells <- lapply(datasets, function(dataset) {
    read_tsv(file.path(raw_dir, dataset, "filtered_feature_bc_matrix/barcodes.tsv.gz"),
             col_names = FALSE) %>% unlist() %>% as.character()
  }) %>% set_names(datasets)
  saveRDS(valid_cells, valid_cell_path)
} else {
  valid_cells <- readRDS(valid_cell_path)
}  

# Create reference seurat objects
refs_path <- file.path(output_dir, "refs.rds")
if (!file.exists(refs_path)) {
  refs <- pbmclapply(datasets, function(dataset) {
    nickname <- str_glue("{dataset}_p1")
    f <- str_glue("{nickname}.csv")
    matrices_to_seurat(file.path(ds_dir, f),
                       nickname = nickname,
                       cells = valid_cells[[dataset]],
                       dataset = dataset,
                       data_type = "reference",
                       condition = "control")
  }, mc.cores = cores) %>% set_names(datasets)
  saveRDS(refs, refs_path)
} else {
  refs <- readRDS(refs_path)
}  


# Compare imputed data with original data
out_all <- list()
out_umap_all <- list()
for (dataset in datasets) {
 
  message("Comparing ", dataset, " >>>>>>>>>>")

  if (dataset == "ds2") {
    combinations <- expand.grid(dataset, ds_cells_props)
  } else if (dataset == "ds3") {
    combinations1 <- expand.grid(dataset, ds_reads_props)
    combinations2 <- expand.grid(dataset, ds_cells_props)
    combinations <- rbind(combinations1, combinations2)
  } else {
    combinations <- expand.grid(dataset, ds_reads_props)
  }
  
  filenames <- str_glue("{combinations$Var1}_{combinations$Var2}")

  querys <-  pbmclapply(filenames, function(f) {
    data_path <- file.path(data_dir, f)
    is_by_cells <-  str_split(basename(f), "_")[[1]][2] %in% ds_cells_props

    new_filename <- str_glue("{team}_{f}.csv")
    if (is_by_cells) { # just take one of replicates if downsampled by cells
      imp <- fread(paste0(data_path, "_n", 1, ".csv")) %>% tibble::column_to_rownames("V1")
      fwrite(imp, new_filename, row.names = TRUE)
    } else {
      new_data <- aggregate_output_by_mean(data_path_list = paste0(data_path, "_n", 1:3, ".csv"),
                                           new_filename = new_filename)
    }
    
    if (file.exists(str_glue("seurat_objects/{team}_{f}.rds"))) {
      so <- readRDS(str_glue("seurat_objects/{team}_{f}.rds"))
    } else {
      so <- matrices_to_seurat(new_filename,
                               nickname = str_glue("{team}_{f}"),
                               cells = valid_cells[[dataset]],
                               dataset = dataset,
                               data_type = "query",
                               condition = "perturb")
    }
    return(so)

  }, mc.cores = length(filenames)) %>% set_names(filenames)
  
  umap_info <- pbmclapply(querys, function(query) {
    # get library size factor
    sf_df <- FetchData(query, vars = "nCount_RNA")
    # create umap info table
    query[["umap"]]@cell.embeddings %>% 
      as.data.frame() %>%
      tibble::rownames_to_column("cell_barcode") %>% 
      mutate(query = query@project.name,
             nCount_RNA = sf_df$nCount_RNA[match(cell_barcode, rownames(sf_df))])
  }, mc.cores = length(filenames)) %>% set_names(filenames)
  out_umap_all[[dataset]] <- umap_info
  
  # out <- pbmclapply(seq_along(querys), function(i) {
  #   compare_so_to_ref(querys[[i]], refs[[dataset]])
  # }, mc.cores = length(filenames)) %>% set_names(filenames)
  # out_all[[dataset]] <- out
}

out_path <-  file.path(output_dir, str_glue("{team}_comparisons.rds"))
saveRDS(out_all, out_path)

out_umap_path <-  file.path(output_dir, str_glue("{team}_umap.rds"))
saveRDS(out_umap_all, out_umap_path)

# if (!file.exists(out_path)) {
#   saveRDS(out_all, out_path)
# } 
# out_all <- readRDS(out_path)


# # Clustering comparison ami vs nmi
# message("Generating clustering comparison plots >>>>>>>>>>")
# clustering_plot_path <- file.path(output_dir, str_glue("{team}_clustering_knn_dropouts_comparison.pdf"))
# pdf(clustering_plot_path)
# for (dataset in datasets) {
#   print(dataset)
#   clust_df <- do.call("rbind",
#                       unlist(lapply(X = out_all[[dataset]],
#                                     FUN =  "[",
#                                     'cluster_metrics'),
#                               recursive = FALSE)) 
  
#   p1 <- ggplot(clust_df, aes(x = ari, y = nmi))+
#     geom_point(aes(size = resolution))+
#     ylim(0,1)+
#     geom_smooth(formula = y ~ x, method = "lm")+
#     theme_bw()+ 
#     facet_wrap(~query_name, scales = "free", ncol = ifelse(dataset == "ds3", 2, 1)) +
#     ggtitle("Leiden clustering overlap") +
#     labs(caption = str_glue("Dataset: {dataset}\n Team: {team}"))
  
#   clust_df2 <- clust_df %>%
#     gather(c(ari, nmi),
#           key = 'metric',
#           value = 'value')

#   p2 <- ggplot(clust_df2, aes(x = resolution, y = value, color = query_name))+
#     geom_line()+
#     geom_point()+
#     facet_wrap(~metric)+
#     theme_bw()+
#     ggtitle("Dot/line plot of clustering ari/nmi") +
#     labs(caption = str_glue("Dataset: {dataset}\n Team: {team}"))
  
#   p3 <- ggplot(clust_df2, aes(x = metric, y = value, fill = query_name))+
#     geom_boxplot()+
#     ylim(0,1) +
#     theme_bw()+
#     ggtitle("Boxplot of clustering ari/nmi") +
#     labs(caption = str_glue("Dataset: {dataset}\n Team: {team}"))
  
#   print(p1)
#   print(p2)
#   print(p3)
  
#   knn_df <- do.call('rbind',
#                     unlist(lapply(X = out_all[[dataset]],
#                                   FUN =  "[",
#                                   'knn_tibble'),
#                             recursive = FALSE))

#   p4 <- ggplot(knn_df, aes(x = query_name, y = jaccard_shared, fill = query_name))+
#     geom_violin(draw_quantiles = c(0.25,0.75))+
#     stat_summary(fun = 'mean',
#                 geom = 'point',
#                 color = 'black')+
#     theme_bw()+
#     ggtitle("Boxplot of KNN jaccard") +
#     labs(caption = str_glue("Dataset: {dataset}\n Team: {team}"))

#   print(p4)

#   dropouts <- do.call('rbind',
#                       unlist(lapply(X = out_all[[dataset]],
#                                     FUN =  "[",
#                                     'drop_outs'),
#                             recursive = FALSE))
#   dropout_rates <- dropouts %>%
#     mutate(n_events = tp+fp+tn+fn,
#           tp_rate = tp/n_events,
#           fp_rate = fp/n_events,
#           tn_rate = tn/n_events,
#           fn_rate = fn/n_events) %>%
#     gather(-c(query_name,
#               reference_name),
#           key = 'metric',
#           value = 'score') %>%
#     filter(!metric %in% c('tp', 'fp', 'tn', 'fn', 'n_events'))
  
#   p5 <- ggplot(dropout_rates, aes(x = metric, y = score, color = query_name))+
#     geom_point(size = 4)+
#     geom_line(aes(group = query_name))+
#     theme_bw()+
#     ggtitle("Dot/line plot for drop out") +
#     labs(caption = str_glue("Dataset: {dataset}\n Team: {team}"))
#   print(p5)
# }
# dev.off()

# Log fold-change comparisons
# message("Generating IFC comparison plots >>>>>>>>>>")
# ifc_plot_path <- file.path(output_dir, str_glue("{team}_log_fold_change_comparison.pdf"))
# pdf(ifc_plot_path)
# groups <- list(c("ds1b", "ds1a"), c("ds1d", "ds1c"))
# lfc_res <- list()
# for (g in groups) {
#   controls <- c(str_glue("{team}_{g[2]}_{ds_reads_props}"), str_glue("{g[2]}_p1"))
#   perturbs <- c(str_glue("{team}_{g[1]}_{ds_reads_props}"), str_glue("{g[1]}_p1"))
#   comparisons_to_make <- data.frame(perturb_nicknames = perturbs,
#                                     control_nicknames = controls)

#   lfcs <- pbmclapply(1:nrow(comparisons_to_make), function(i) {
#     lfc_comparison(comparisons_to_make[i, ])
#   }, mc.cores = cores) %>% set_names(gsub(str_glue("{team}_"), "", comparisons_to_make$control_nicknames))
  
  
#   lfc_10k <- full_join(x = lfcs[[1]],
#                         y = lfcs[[4]],
#                         by = "gene",
#                         suffix = c("_10k", "_1")) %>%
#     mutate(avg_log2FC_1 = ifelse(is.na(avg_log2FC_1), 0, avg_log2FC_1),
#            avg_log2FC_10k = ifelse(is.na(avg_log2FC_10k), 0, avg_log2FC_10k),
#            team = team, comparsion_group = paste(g, collapse = " vs "))

#   lfc_20k <- full_join(x = lfcs[[2]],
#                       y = lfcs[[4]],
#                       by = "gene",
#                       suffix = c("_20k", "_1")) %>%
#     mutate(avg_log2FC_1 = ifelse(is.na(avg_log2FC_1), 0, avg_log2FC_1),
#            avg_log2FC_20k = ifelse(is.na(avg_log2FC_20k), 0, avg_log2FC_20k),
#            team = team, comparsion_group = paste(g, collapse = " vs "))

#   lfc_50k <- full_join(x = lfcs[[3]],
#                       y = lfcs[[4]],
#                       by = "gene",
#                       suffix = c("_50k", "_1")) %>%
#     mutate(avg_log2FC_1 = ifelse(is.na(avg_log2FC_1), 0, avg_log2FC_1),
#            avg_log2FC_50k = ifelse(is.na(avg_log2FC_50k), 0, avg_log2FC_50k),
#            team = team, comparsion_group = paste(g, collapse = " vs "))

#   g_lfc_res <- list(lfc_10k = lfc_10k, lfc_20k = lfc_20k, lfc_50k = lfc_50k)
#   lfc_res <- append(lfc_res, g_lfc_res)

  # p1 <- ggplot(lfc_10k, aes(x = avg_log2FC_1, y = avg_log2FC_10k))+
  #   geom_point()+
  #   geom_smooth(formula = y ~ x, method = 'lm')+
  #   theme_bw()+
  #   coord_equal()+
  #   labs(title = 'LFC correlation: 100% vs 10k reads',
  #       caption = str_glue("Perturb ({g[2]}) vs Control ({g[1]}) \n Team: {team}"),
  #       subtitle = paste0('Pearson cor: ', cor(lfc_10k$avg_log2FC_1, lfc_10k$avg_log2FC_10k, use = 'complete.obs')))

  # p2 <- ggplot(lfc_20k, aes(x = avg_log2FC_1, y = avg_log2FC_20k))+
  #   geom_point()+
  #   geom_smooth(formula = y ~ x, method = 'lm')+
  #   theme_bw()+
  #   coord_equal()+
  #   labs(title = 'LFC correlation: 100% vs 20k reads',
  #       caption = str_glue("Perturb ({g[2]}) vs Control ({g[1]}) \n Team: {team}"),
  #       subtitle = paste0('Pearson cor: ', cor(lfc_20k$avg_log2FC_1, lfc_20k$avg_log2FC_20k, use = 'complete.obs')))

  # p3 <- ggplot(lfc_50k, aes(x = avg_log2FC_1, y = avg_log2FC_50k))+
  #   geom_point()+
  #   geom_smooth(formula = y ~ x, method = 'lm')+
  #   theme_bw()+
  #   coord_equal()+
  #   labs(title = 'LFC correlation: 100% vs 50k reads',
  #       caption = str_glue("Perturb ({g[2]}) vs Control ({g[1]}) \n Team: {team}"),
  #       subtitle = paste0('Pearson cor: ', cor(lfc_50k$avg_log2FC_1, lfc_50k$avg_log2FC_50k, use = 'complete.obs')))
  
  # print(p1)
  # print(p2)
  # print(p3)
# }
# dev.off()

# saveRDS(lfc_res, file.path(output_dir, str_glue("{team}_lfc_results.rds")))

# upload to synapse
# message("Uploading plots to synpase >>>>>>>>>>")
# system(str_glue("synapse store {clustering_plot_path} --parentId syn52812201"))
# system(str_glue("synapse store {ifc_plot_path} --parentId syn52812201"))

# Remove temp files
message("Removing temp files >>>>>>>>>>")
rm_files <- list.files(".", pattern = str_glue("{team}_*"), full.names = TRUE)
unlink(rm_files)
# rm_files2 <- list.files("seurat_objects", pattern = str_glue("{team}_*"), full.names = TRUE)
# unlink(rm_files2)
gc();gc()
