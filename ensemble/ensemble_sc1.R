suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(data.table)
  library(patchwork)
  library(ggsci)
})

# added in .Rprofile
# source("utils.R")
# source("bootstrap_funcs.R")
# source("plot_funcs.R")
# source("metrics.R")

ncores <- parallel::detectCores() - 1


# Set up synapse ----------------------------------------------------------
reticulate::use_condaenv('synapse')
synapseclient <- reticulate::import('synapseclient')
syn <- synapseclient$Synapse()
syn$login(silent = TRUE)


# Download submissions ----------------------------------------------------
task_n <- "task1"
view_id <- "syn51157023"
eval_id <- "9615023"
gs_id <- "syn34612394"
metrics_lookup <- c("nrmse_score", "spearman_score")

# query submissions
all_sub_df <- get_ranked_submissions(syn, view_id, eval_id, "private")

# label model names
baseline_magic <- "9732066"
baseline_deepimpute <- "9732074"
top_performer <- all_sub_df$id[1]
all_sub_df <- all_sub_df %>% mutate(model_name = case_when(id == baseline_magic ~ "Baseline MAGIC",
                                                           id == baseline_deepimpute ~ "Baseline DeepImpute",
                                                           TRUE ~ as.character(team)))


# Ensemble analysis ----------------------------------------------------

# exclude the baseline to run ensemble
sub_df <- all_sub_df %>% filter(!id %in% c(baseline_magic, baseline_deepimpute))

# download the ground truth
gs_path <- syn$get(gs_id)["path"]
all_gs <- readRDS(gs_path)

# prepare
all_scores <- data.table()
n <- 0
temp_dir <- "temp"
dir.create(temp_dir, showWarnings = FALSE)

message("!!!!!!!!!!!!! Ensemble steps below requires high I/O")
message("!!!!!!!!!!!!! and high computing due to the large number of scrna data")
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
  pred_dir <- file.path(temp_dir, str_glue("{team}_{sub_id}"))
  if (!dir.exists(pred_dir)) {
    pred_path <- syn$get(pred_id)$path
    dir.create(pred_dir, showWarnings = FALSE, recursive = TRUE)
    # decompress the prediction files
    untar(pred_path, exdir = pred_dir)
  }
  # get all prediction file names
  pred_files <- list.files(file.path(pred_dir, "output"), pattern = "_imputed.csv")
  
  message("Calculating scores on all prediction files ...")
  if (n > 0) message("Ensembling by mean ...")
  
  scores_df <- pbmcapply::pbmclapply(pred_files, function(pred_file) {
    
    # detect file prefix used to read gs
    info <- strsplit(pred_file, "_")[[1]]
    prefix <- info[1]
    prop <- info[2]
    
    # read pred file and convert to long format
    pred_path <- file.path(pred_dir, "output", pred_file)
    pred_data <- fread(pred_path, verbose = FALSE) %>%
      melt.data.table(id.vars = 1) %>% 
      setnames(c("genes", "cells", "values"))
    
    # combine by taking mean of values
    if (n > 0) {
      ref_data <- readRDS(file.path(temp_dir, str_glue("ref_{pred_file}.rds")))
      pred_data <- rbindlist(list(ref_data, pred_data), use.names = TRUE) %>%
        .[, .(values = mean(values)), by = .(cells, genes)]
    }
    
    # save new reference data
    saveRDS(pred_data, file.path(temp_dir, str_glue("ref_{pred_file}.rds")))
    
    # convert back to gene expr matrix
    pred_data <- pred_data %>%
      dcast(genes ~ cells, value.var = "values") %>%
      tibble::column_to_rownames("genes")
    
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
    return(
      data.table(dataset = pred_file,
                 nrmse_score = nrmse_score,
                 spearman_score = abs(spearman_score),
                 ensemble_label = label)
    )
  }, mc.cores = ncores, ignore.interactive = TRUE) %>%
    rbindlist(use.names = TRUE)
  
  if (nrow(scores_df) != length(pred_files)) stop("Ensemble process failed on ", sub_id, ".\n",
                                                  "If still not work, try to do ensemble one submission at a time manaully.")
  all_scores <- rbindlist(list(all_scores, scores_df), use.names = TRUE)
  
  # save all scores in case partially failed
  # could add cache step, but no time for now
  saveRDS(all_scores, "ensemble_all_scores_sc1.rds")
  n <- n + 1
  gc(); gc()
}

# remove the temp dir after the ensemble process completes
# unlink(temp_dir, recursive = TRUE, force = TRUE)

# read the ensembled scores
all_scores <- readRDS("ensemble_all_scores_sc1.rds")

# add back baseline scores
baseline_scores <- all_sub_df %>% 
  filter(id %in% c(baseline_magic, baseline_deepimpute)) %>%
  get_scores(syn, sub_df = .) %>%
  mutate(ensemble_label = all_sub_df$model_name[match(id, all_sub_df$id)]) %>%
  select(dataset, nrmse_score, spearman_score, ensemble_label) 

all_scores <- bind_rows(all_scores, baseline_scores)


# Rank new submissions --------------------------------------------------------
# rank submissions all ensemble models and baselines
# will be used to plot later
rank_df <- all_scores %>%
  rank_submissions(metrics_lookup[1], metrics_lookup[2], "ensemble_label") %>%
  mutate(team = gsub("\\d? \\+ (.*)$", "\\1", ensemble_label),
         ensemble_label = factor(ensemble_label, levels = c("0 + GOAL_LAB", 
                                                            "1 + DLS5", 
                                                            "2 + LDExplore",
                                                            "Baseline MAGIC",
                                                            "3 + Metformin-121",
                                                            "Baseline DeepImpute")))


# Plotting ranks over ensembled models ------------------------------------
my_labels <- c(primary_rank = "NRMSE", secondary_rank = "Spearman Correlation")
my_colors <- c(primary_rank = "#709AE1E5", secondary_rank = "#FED439E5")
rank_line_p <- rank_df %>% 
  gather("metrics", "ranks", c(2:3)) %>%
  mutate(ranks = 1 / ranks) %>%
  ensemble_ranks_line(ensemble_label, ranks, metrics) +
  scale_color_manual(values = my_colors, labels = my_labels) +
  labs(x = NULL, y = "1 / Ranks", color = "Metrics") +
  theme(legend.position = c(0.95, 0.95), legend.justification = c("right", "top"),
        panel.grid.minor.y = element_blank(),
        text = element_text(size = 18),
        axis.title = element_text(size = 20))


# Bootstrapping -----------------------------------------------------------
# bootstrapping the rankings
boot_df <- bootstrap(.data = all_scores,
                     seq_size = length(unique(all_scores$dataset)),
                     .by = "ensemble_label",
                     n_iterations = 1000,
                     seed = 165136,
                     ncores = ncores)

# rank submissions for all bootstraps
boot_rank_df <- boot_df %>%
  nest(scores = c(metrics_lookup[1], metrics_lookup[2], dataset, ensemble_label), .by = bs_n) %>%
  mutate(ranks = parallel::mclapply(scores, rank_submissions, metrics_lookup[1], metrics_lookup[2], "ensemble_label", mc.cores = ncores)) %>%
  unnest(cols = ranks) %>%
  select(-scores)

# compute bayes factor of each new model against top performer
ref_ranks <- boot_rank_df %>% filter(ensemble_label == "0 + GOAL_LAB")
bf_df <- boot_rank_df %>%
  group_by(ensemble_label) %>%
  mutate(primary_bf = bayes_factor(primary_rank, ref_ranks$primary_rank),
         secondary_bf = bayes_factor(secondary_rank, ref_ranks$secondary_rank)) %>%
  gather("metrics", "ranks", c(primary_rank, secondary_rank)) %>%
  mutate(ensemble_label = factor(ensemble_label, 
                                 levels = c("0 + GOAL_LAB", 
                                            "1 + DLS5", 
                                            "2 + LDExplore",
                                            "Baseline MAGIC",
                                            "3 + Metformin-121",
                                            "Baseline DeepImpute")))



# Plotting ----------------------------------------------------------------
# against top performer
p_top1 <- bf_df %>%
  filter(metrics == "primary_rank") %>%
  mutate(ranls = 1 / ranks) %>%
  bootstrap_boxplot(ensemble_label, ranks, primary_bf, 
                    bf_cutoffs=c(5, 30), ref_label = "Ref: Top Performer") +
  labs(x = NULL, y = "1 / (Bootstrapped Ranks of NRMSE)", color = "Bayes Factor")

p_top2 <- bf_df %>%
  filter(metrics == "secondary_rank") %>%
  mutate(ranls = 1 / ranks) %>%
  bootstrap_boxplot(ensemble_label, ranks, primary_bf, 
                    bf_cutoffs=c(5, 30), ref_label = "Ref: Top Performer") +
  labs(x = NULL, y = "1 / (Bootstrapped Ranks of Spearman Correlation)", color = "Bayes Factor")

# saving all plots
pdf(file="sc1_ensemble_analysis.pdf", width = 16, height = 10)
rank_line_p;
p_top1 / p_top2 + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "top", legend.direction = "horizontal")
dev.off()