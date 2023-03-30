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
task_n <- "task2"
view_id <- "syn51157023"
eval_id <- "9615024"
gs_id <- "syn35294386"
all_sub_df <- get_ranked_submissions(syn, view_id, eval_id, "private")

baseline_macs2 <- "9732044"
top_performer <- all_sub_df$id[1]
all_sub_df <- all_sub_df %>% mutate(model_name = case_when(id == baseline_macs2 ~ "Baseline Macs2",
                                                           TRUE ~ as.character(team)))


# exclude the baseline for ensemble analysis
sub_df <- all_sub_df %>% filter(!id %in% c(baseline_macs2))

# download the ground truth
gs_path <- syn$get(gs_id)["path"]
all_gs <- readRDS(gs_path)

# iterate each submission to perform bootstrapping
base_data_list <- list()
all_scores <- tibble()
n <- 0

while (n < nrow(sub_df)) { 
  
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
  temp_dir <- str_glue("sc2_{label}")
	dir.create(temp_dir, showWarnings = FALSE)

	pred_path <- syn$get(pred_id)$path
	pred_data <- fread(pred_path)

  message("Calculating scores on all prediction files (randomly pick to reduce memory usuage)...")
  if (n > 0) message("Ensembling by max ...")
  
  scores_res <- pbmcapply::pbmclapply(, function(i) {
    
    pred_file <- pred_files[i]
    # detect file prefix used to read gs
    info <- strsplit(pred_file, "\\.")[[1]]
    prefix <- info[1]

    # read gs
    gs <- all_gs$gs_data[[sub_phase]][[prefix]]
    gs_ranked_filtered <- gs[["gs_ranked_filtered"]]
    gs_sort <- gs[["gs_sort"]]
    gs_sort_ubi <- gs[["gs_sort_ubi"]]
    gs_sort_tss <- gs[["gs_sort_tss"]]
    gs_sort_csp <- gs[["gs_sort_csp"]]
  
    # default scores
    jaccard_similarity <- 0
    recall_ubiquitous <- 0
    recall_tss <- 0
    recall_cellspecific <- 0
    summed_score <- 0

    # read prediction
    # pred_data <- data.table::fread(file.path(temp_dir, "output", pred_file), data.table = FALSE, verbose = FALSE)
    
    if (n > 0) {
      df0 <- base_data_list[[pred_file]]
      df1 <- pred_data
      merged_df <- rbind(df0, df1)

      # update pred_data with combined values - much faster using data.table
      setDT(merged_df)
      # combine by taking max of values
      avg_df <- merged_df[, .(avgs = max(values)), by = .(chr, start, end)]
      pred_df <- avg_df %>%
        spread(cells, avgs) %>%
        tibble::column_to_rownames("genes")
    }
    fst::write_fst(pred_df, file.path(temp_dir, pred_file))
    if (nrow(pred_data) > 0) {
      
      pred_data <- pred_data %>% filter(grepl("chr(\\d|X|Y)", chr)) %>% select(1:3) %>% as.data.frame() %>% setNames(c("chr", "start", "end"))
      pred_sort <- bedr.merge.region(bedr.sort.region(pred_data, check.zero.based = FALSE, check.chr = FALSE, check.valid = FALSE, verbose = FALSE), verbose = FALSE)

      # caculate scores
      j_result <- jaccard(gs_sort, pred_sort, check.merge = TRUE, check.chr = FALSE, check.sort = FALSE, check.valid = FALSE, verbose = FALSE)
      jaccard_similarity <- as.numeric(j_result$jaccard)
      recall_ubiquitous <- category_recall(gs_sort_ubi, pred_sort)
      recall_tss <- category_recall(gs_sort_tss, pred_sort)
      recall_cellspecific <- category_recall(gs_sort_csp, pred_sort)

      summed_score <- jaccard_similarity + recall_ubiquitous + recall_tss

      # report cellspecific for ds2, but not added to summed score
      if (prefix == "ds1") summed_score <- summed_score + recall_cellspecific
    }

    # collect scores for each test case
    return(
      tibble(
        dataset = pred_file,
        jaccard_similarity = jaccard_similarity,
        recall_ubiquitous = recall_ubiquitous,
        recall_tss = recall_tss,
        recall_cellspecific = recall_cellspecific,
        summed_score = summed_score
      )
    )
    return(list(data = pred_data, scores = scores_df))
  }, mc.cores = ncores, ignore.interactive = TRUE)

  scores_df <- lapply(scores_res, `[[`, 2) %>% bind_rows()
  stopifnot(nrow(scores_df) == length(pred_files))
  all_scores <- rbind(all_scores, scores_df)
  fst::write_fst(all_scores, "ensemble_all_scores_sc2.fst")

  # update base data after ensembled with mean for next submission
  base_data_list <- lapply(scores_res, `[[`, 1) %>% setNames(pred_files)
  n <- n + 1
  gc(); gc()
  # unlink(temp_dir, recursive = TRUE)
}



