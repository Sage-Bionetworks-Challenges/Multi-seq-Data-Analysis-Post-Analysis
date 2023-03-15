library(ggplot2)
library(dplyr)
library(stringr)

source("metrics.R")
source("utils.R")
ncores <- parallel::detectCores() - 1

# download the submissions
task_n <- "task1"
view_id <- "syn51157023"
gs_id <- "syn34612394"
# gs_ids <- list(task1 = "syn34612394", task2 = "syn35294386")
sub_df <- get_submissions(view_id, gs_id)

# download the ground truth
gs_path <- syn$get(task_gs_id)["path"]
all_gs <- readRDS(gs_path)


task_subs <- sub_df[sub_df$task == task_n,]

# iterate each submission to perform bootstrapping
all_scores <- tibble()
for (sub_n in 1:nrow(task_subs)) { # use for loop to save computation for bootstrapping
  
  task_pred_id <- task_subs$prediction_fileid[sub_n]
  task_team <- task_subs$team[sub_n]
  # download the submission's prediction file
  # it will take some time, since each prediction file is a large tarball
  pred_path <- syn$get(task_pred_id)$path
  
  # create a temp dir to store all prediction files for one submission
  pred_out_dir <- str_glue("{task_team}_{task_subs$id[i]}")
  dir.create(pred_dir, showWarnings = FALSE)
  
  # decompress the prediction files into the temp dir
  message("Untaring submission to ", pred_out_dir, " ...")
  untar(pred_path, exdir = pred_out_dir)
  
  # bootstrap 1000 times
  bs_n <- 1000
  pred_files <- list.files(pred_dir, pattern = "_imputed.csv")
  bs_indices <- bootstrap(seq_size = length(pred_files), n_iterations = bs_n, seed = 98109)
  
  # collect all bootstrapping scores for each submission
  bs_scores <- parallel::mclapply(sequence(bs_n), function(n) {
    bs_indice <- bs_indices[[n]]
    
    scores_df <- lapply(pred_files[bs_indice], function(pred_file) {
      
      # detect file prefix used to read gs
      info <- strsplit(pred_file, "_")[[1]]
      prefix <- info[1]
      prop <- info[2]
      
      # read prediction
      pred_path <- file.path(pred_dir, pred_file)
      pred_data <- data.table::fread(pred_path, data.table = FALSE, verbose = FALSE) %>% tibble::column_to_rownames("V1")
      pred_data <- Seurat::NormalizeData(pred_data, verbose = FALSE)
      
      # read gs
      gs_data <- all_gs$gs_data[[prefix]]
      
      use_pseudobulk <- prop %in% c("p00625", "p0125", "p025")
      
      # prepare the scrna data for evaluation
      eval_data <- .prepare(
        true = gs_data,
        pred = pred_data,
        pseudobulk = use_pseudobulk,
        proportion = prop
      )
      
      # scoring
      nrmse_score <- calculate_nrmse(eval_data, pseudobulk = use_pseudobulk)
      spearman_score <- calculate_spearman(eval_data, pseudobulk = use_pseudobulk)
      
      # collect scores for each test case
      return(
        tibble(
          dataset = pred_file,
          nrmse_score = nrmse_score,
          spearman_score = abs(spearman_score),
          bs_n = n
        )
      )
    }) %>% bind_rows()
    
    # correct direction of nrmse scores for ranking
    bs_scores$nrmse_score <- -bs_scores$nrmse_score
    
    # rank all submissions
    rank_df <- rank_submission(scores = socres_df, 
                               dataset_col = "dataset",
                               bs_col = "bs_n",
                               primary = "nrmse_score", 
                               secondary = "spearman_score")
    
    # remove scored submission prediction files and temp dir to clean space
    unlink(pred_out_dir, recursive = TRUE) 
  }, mc.cores = ncores) %>% bind_rows()

  
  # append the each submission's bootstrapping results
  all_scores <- rbind(all_scores, bs_scores)
}

# tbd for bayes factor



