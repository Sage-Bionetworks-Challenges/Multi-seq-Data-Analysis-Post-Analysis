# download the submissions
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

source("metrics.R")
source("utils.R")
ncores <- parallel::detectCores() - 1

# set up synapse
reticulate::use_condaenv('synapse')
synapseclient <- reticulate::import('synapseclient')
syn <- synapseclient$Synapse()
syn$login(silent = TRUE)


# download the submissions
task_n <- "task1"
view_id <- "syn51157023"
eval_id <- "9615023"
gs_id <- "syn34612394"
# gs_ids <- list(task1 = "syn34612394", task2 = "syn35294386")
# each team should be ranked only on the best submission, aka one submission per team from this table
sub_df <- get_ranked_submissions(syn, view_id, eval_id, "private")

# download the ground truth
gs_path <- syn$get(gs_id)["path"]
all_gs <- readRDS(gs_path)

# iterate each submission to perform bootstrapping
all_scores <- tibble()

for (sub_n in 1:nrow(sub_df)) { # use for loop to save computation for bootstrapping
  
  sub_id <- sub_df$id[sub_n]
  pred_id <- sub_df$prediction_fileid[sub_n]
  team <- sub_df$team[sub_n]

  message(str_glue("Scoring {team} ..."))
  # download the submission's prediction file
  # it will take some time, since each prediction file is a large tarball
  
  temp_dir <- str_glue("{team}_{sub_id}")

  if (!file.exists(temp_dir)) {
    pred_path <- syn$get(pred_id)$path
    # create a temp dir to store all prediction files for one submission
    dir.create(temp_dir, showWarnings = FALSE)
    
    # decompress the prediction files into the temp dir
    message("Untaring submission to ", temp_dir, " ...")
    untar(pred_path, exdir = temp_dir)
  }
  pred_files <- list.files(file.path(temp_dir, "output"), pattern = "_imputed.csv")

  # collect all bootstrapping scores for each submission
  scores_df <- lapply(seq_along(pred_files), function(i) {
    
    pred_file <- pred_files[i]
    # detect file prefix used to read gs
    info <- strsplit(pred_file, "_")[[1]]
    prefix <- info[1]
    prop <- info[2]
      
    # read prediction
    pred_path <- file.path(temp_dir, "output", pred_file)
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

    # bootstrap 1000 times
    # make indices of each file are the same for all submissions
    bs_n <- 10
    bs_indices <- boot_indices(seq_size = nrow(eval_data$pred), n_iterations = bs_n, seed = 98109)

    bs_scores <- pbmcapply::pbmclapply(seq_along(bs_indices), function(n) {
      
      bs_indice <- bs_indices[[n]]
      eval_data$true <- eval_data$true[bs_indice, ]
      eval_data$pred <- eval_data$pred[bs_indice, ]
      
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
    }, mc.cores = ncores) %>% bind_rows()
    
    # rank_df <- boot_rank_submission(scores = socres_df, 
    #                                 dataset_col = "dataset",
    #                                 bs_col = "bs_n",
    #                                 primary = "nrmse_score", 
    #                                 secondary = "spearman_score")
    print(str_glue("Done {i} / {length(pred_files)}"))

    return(bs_scores)
  }) %>% 
    bind_rows() %>%
    mutate(id = sub_id, team = team)
  
  # remove scored submission prediction files and temp dir to clean space
  # unlink(temp_dir, recursive = TRUE)
  # append the each submission's bootstrapping results
  all_scores <- rbind(all_scores, scores_df)
}

saveRDS(all_scores, "all_scores.rds")

# rank submissions for all bootstraps
rank_df <- boot_df %>%
  nest(scores = c(metrics_lookup[1], metrics_lookup[2], dataset, id, team), .by = bs_n) %>% 
  mutate(ranks = parallel::mclapply(scores, rank_submissions, metrics_lookup[1], metrics_lookup[2], mc.cores = ncores)) %>%
  unnest(cols = ranks) %>%
  select(-scores) %>%
  mutate(model_name = sub_df$model_name[match(id, sub_df$id)])

# compute bayes factor between each model
ref_ids <- c(baseline_magic, baseline_deepimpute, top_performer)
ref_names <- c("baseline_magic", "baseline_deepimpute", "top_performer")
bf_df <- lapply(seq_along(ref_ids), function(i) {
  ref_ranks <- rank_df %>% filter(id == ref_ids[i])
  rank_df %>%
    group_by(id) %>%
    mutate(primary_bf = bayes_factor(primary_rank, ref_ranks$primary_rank),
           secondary_bf = bayes_factor(secondary_rank, ref_ranks$secondary_rank),
           ref_model = ref_names[i])
  }) %>% 
  bind_rows() %>%
  gather("metrics", "ranks", c(primary_rank, secondary_rank))