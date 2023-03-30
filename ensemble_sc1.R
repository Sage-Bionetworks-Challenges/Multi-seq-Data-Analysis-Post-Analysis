# download the submissions
suppressPackageStartupMessages({
	library(ggplot2)
	library(dplyr)
	library(tidyr)
	library(stringr)
	library(data.table)
})

source("metrics.R")
source("utils.R")
source("score_scrna.R")
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
all_sub_df <- get_ranked_submissions(syn, view_id, eval_id, "private")

baseline_magic <- "9732066"
baseline_deepimpute <- "9732074"
top_performer <- all_sub_df$id[1]
all_sub_df <- all_sub_df %>% mutate(model_name = case_when(id == baseline_magic ~ "Baseline MAGIC",
																													 id == baseline_deepimpute ~ "Baseline DeepImpute",
																													 TRUE ~ as.character(team)))

# exclude the baseline for ensemble analysis
sub_df <- all_sub_df %>% filter(!id %in% c(baseline_magic, baseline_deepimpute))

# download the ground truth
gs_path <- syn$get(gs_id)["path"]
all_gs <- readRDS(gs_path)

# iterate each submission to perform bootstrapping
base_data_list <- list()
all_scores <- tibble()
n <- 0

while(n < nrow(sub_df)) { 
  
  sub_n <- n + 1
  sub_id <- sub_df$id[sub_n]
  pred_id <- sub_df$prediction_fileid[sub_n]
  team <- as.character(sub_df$team[sub_n])
  label <- str_glue("{n} + {team}")
	message("\n")
	message("------------------------------------------")
	message("Add predictions from ", team)
	message("------------------------------------------")
  
	message("Retrieving prediction files ...")
	# create a temp dir to store all prediction files for one submission
	temp_dir <- str_glue("{team}_{sub_id}")
	if (!file.exists(temp_dir)) {
		pred_path <- syn$get(pred_id)$path
		dir.create(temp_dir, showWarnings = FALSE)
		# decompress the prediction files
		untar(pred_path, exdir = temp_dir)
	}
	# get all prediction file names
  pred_files <- list.files(file.path(temp_dir, "output"), pattern = "_imputed.csv")

  message("Calculating scores on all prediction files (randomly pick to reduce memory usuage)...")
  if (n > 0) message("Ensembling by mean ...")

  set.seed(1234)
	indice_list <- split(
		sample(seq_along(pred_files)),
		cut(seq_along(pred_files), 5, labels = FALSE)
	)

  scores_res <- list()
	for (indices in indice_list) {
    res <- pbmcapply::pbmclapply(indices, function(i) {
      
      pred_file <- pred_files[i]
      # detect file prefix used to read gs
      info <- strsplit(pred_file, "_")[[1]]
      prefix <- info[1]
      prop <- info[2]
        
      pred_path <- file.path(temp_dir, "output", pred_file)
      pred_data <- fread(pred_path, data.table = FALSE, verbose = FALSE) %>% tibble::column_to_rownames("V1")

      if (n > 0) {
        # combine by taking mean of values
        df0 <- base_data_list[[pred_file]] %>% tibble::rownames_to_column("genes") %>% gather("cells", "values", -genes)
        df1 <- pred_data %>% tibble::rownames_to_column("genes") %>% gather("cells", "values", -genes)
        merged_df <- rbind(df0, df1) 

        # update pred_data with combined values - much faster using data.table
        setDT(merged_df)
        avg_df <- merged_df[, .(avgs = mean(values)), by = .(cells, genes)]
        pred_df <- avg_df %>%
          spread(cells, avgs) %>%
          tibble::column_to_rownames("genes")
      }

      # normalize data
      norm_pred_data <- Seurat::NormalizeData(pred_data, verbose = FALSE)

      # read gs
      gs_data <- all_gs$gs_data[[prefix]]

      use_pseudobulk <- prop %in% c("p00625", "p0125", "p025")

      # prepare the scrna data for evaluation
      eval_data <- .prepare(
        true = gs_data,
        pred = norm_pred_data,
        pseudobulk = use_pseudobulk,
        proportion = prop
      )

      # scoring
      nrmse_score <- calculate_nrmse(eval_data, pseudobulk = use_pseudobulk)
      spearman_score <- calculate_spearman(eval_data, pseudobulk = use_pseudobulk)
      
      # collect scores for each test case
      scores_df <- tibble(dataset = pred_file,
                          nrmse_score = nrmse_score,
                          spearman_score = abs(spearman_score),
                          ensemble_label = label)
      
      return(list(data = pred_data, scores = scores_df))
    }, mc.cores = ncores, ignore.interactive = TRUE)

    scores_res <- append(scores_res, res)
  }

  scores_df <- lapply(scores_res, `[[`, 2) %>% bind_rows()
  stopifnot(nrow(scores_df) == length(pred_files))
  all_scores <- rbind(all_scores, scores_df)
  fst::write_fst(all_scores, "ensemble_all_scores_sc1.fst")

  # update base data after ensembled with mean for next submission
  base_data_list <- lapply(scores_res, `[[`, 1) %>% setNames(pred_files)
  n <- n + 1
  gc(); gc()
  # unlink(temp_dir, recursive = TRUE)
}



