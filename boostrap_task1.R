library(ggplot2)
library(dplyr)
library(stringr)

source("metrics.R")
source("utils.R")
ncores <- parallel::detectCores() - 1

# download the submissions
view_id <- "syn51157023"
sub_df <- get_submissions(view_id)


for (task in c("task1", "task2")) {
  
  # get all scores
  scores_df <- get_scores(sub_df[sub_df$task == task, ])
  
  # get corresponding metric names and for each task
  if (task == "task1") {
    metrics_lookup <- c("nrmse_score", "spearman_score")
    # correct the direction of nrmse
    scores_df$nrmse_score <- -scores_df$nrmse_score
  } else {
    metrics_lookup <- c("summed_score", "jaccard_similarity")
  }

  # bootstrapping the rankings
  boot_res <- bootstrap(.data = scores_df,
                        seq_size = length(unique(scores_df$dataset)),
                        func = rank_submissions,
                        primary_metric = metrics_lookup[1],
                        secondary_metric = metrics_lookup[2],
                        .by = "id",
                        n_iterations=1000,
                        seed=98109,
                        ncores=ncores) 
}

# tbd for bayes factor



