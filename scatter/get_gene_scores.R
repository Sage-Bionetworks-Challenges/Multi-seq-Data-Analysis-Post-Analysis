source("utils/setup.R")
source("utils/bootstrap_funcs.R")
source("utils/plot_funcs.R")

task_n <- 1
metrics <- metrics_lookup[[task_n]]


# Reading submission data -------------------------------------------------
sub_data <- file.path(data_dir, str_glue("final_submissions_task{task_n}.rds"))
if (!file.exists(sub_data)) source("submission/get_submissions.R")

sub_df <- readRDS(sub_data)


# label model names -------------------------------------------------------
baseline_magic <- "9732066"
baseline_deepimpute <- "9732074"
top_performer <- sub_df$id[1]
sub_df <- sub_df %>% mutate(model_name = case_when(id == baseline_magic ~ "Baseline MAGIC",
                                                   id == baseline_deepimpute ~ "Baseline DeepImpute",
                                                   TRUE ~ as.character(team)))


# Get scores on all genes ----------------------------------------------------------

# download the ground truth
gs_path <- syn$get(gs_id[[task_n]])["path"]
all_gs <- readRDS(gs_path)
  
# prepare
all_scores <- data.table()
temp_dir <- "temp"
dir.create(temp_dir, showWarnings = FALSE)

message("!!!!!!!!!!!!! Ensemble steps below requires high I/O")
message("!!!!!!!!!!!!! and high computing due to the large number of scrna data")
for (n in 1:nrow(sub_df)) { 
  
  sub_n <- n
  sub_id <- sub_df$id[sub_n]
  pred_id <- sub_df$prediction_fileid[sub_n]
  team <- as.character(sub_df$team[sub_n])

  message("\n")
  message("------------------------------------------")
  message("Scoring submissions from ", team)
  message("------------------------------------------")
  
  message("Retrieving prediction files ...")
  pred_dir <- file.path(temp_dir, str_glue("{team}_{sub_id}_task{task_n}"))
  if (!dir.exists(pred_dir)) {
    source(str_glue("submission/get_predictions_task{task_n}.R"))
  }
  
  # get all prediction file names
  pred_files <- list.files(file.path(pred_dir, "output"), pattern = "_imputed.csv")
  
  message("Calculating scores on each gene...")
  scores_df <- pbmcapply::pbmclapply(pred_files, function(pred_file) {
    
    # detect file prefix used to read gs
    info <- strsplit(pred_file, "_")[[1]]
    prefix <- info[1]
    prop <- info[2]
    
    # read pred file and convert to long format
    pred_path <- file.path(pred_dir, "output", pred_file)
    pred_data <- fread(pred_path, verbose = FALSE) %>% tibble::column_to_rownames("V1")
    pred_data <- Seurat::NormalizeData(pred_data, verbose = FALSE)
    
    # read gs
    gs_data <- all_gs$gs_data[[prefix]]
    gs_data <- Seurat::NormalizeData(gs_data, verbose = FALSE)
    
    # prepare the scrna data for evaluation
    use_pseudobulk <- prop %in% c("p00625", "p0125", "p025")
    eval_data <- .prepare(
      true = gs_data,
      pred = pred_data,
      pseudobulk = use_pseudobulk,
      proportion = prop
    )
    
    # return the vector of scores on all genes
    nrmse_genes_score <- calculate_nrmse(eval_data, pseudobulk = use_pseudobulk)
    spearman_genes_score <- calculate_spearman(eval_data, pseudobulk = use_pseudobulk, ncores = ncores / 2)
    
    # collect scores of all genes for each test case
    return(
      data.table(dataset = pred_file,
                 nrmse_score = nrmse_genes_score,
                 spearman_score = abs(spearman_genes_score),
                 genes = rownames(eval_data$true))
    )
  }, mc.cores = ncores, ignore.interactive = TRUE) %>%
    rbindlist(use.names = TRUE) %>%
    mutate(id = sub_id, team = team)
  
  
  all_scores <- rbindlist(list(all_scores, scores_df), use.names = TRUE)
  
  gc(); gc()
}

# will take a while due to large numebr of rows
all_scores <- all_scores %>%
  .[,  c("prefix", "prop", "replicate", "dummy") := tstrsplit(dataset, "_")] %>%
  mutate(model_name = sub_df$model_name[match(all_scores$id, sub_df$id)]) %>%
  mutate(downsampled = if_else(prop %in% c("p00625", "p0125", "p025"), "by_cells", "by_reads"),
         model_name = factor(model_name, levels = unique(model_name)),
         prefix = toupper(prefix),
         prop = factor(case_when(prop == "p00625" ~ "6.25%",
                                 prop == "p0125" ~ "12.5%",
                                 prop == "p025" ~ "25%",
                                 TRUE ~ gsub("p(\\d+k)", "\\1", prop)),
                       levels = c("6.25%", "12.5%", "25%", "10k", "20k", "50k"))) %>%
  select(-dummy)

saveRDS(all_scores, file.path(data_dir, "all_genes_scores_task1.rds"))

