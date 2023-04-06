# download the submissions
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(data.table)
  library(patchwork)
  library(ggsci)
})

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
all_scores <- data.table()
n <- 0
temp_dir <- "temp"
dir.create(temp_dir, showWarnings = FALSE)

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
  
  all_scores <- rbindlist(list(all_scores, scores_df), use.names = TRUE)
  
  n <- n + 1
  gc(); gc()
}
# unlink(temp_dir, recursive = TRUE, force = TRUE)

saveRDS(all_scores, "ensemble_all_scores_sc1.rds")
all_scores <- readRDS("ensemble_all_scores_sc1.rds")

# now add baseline scores
baseline_scores <- all_sub_df %>% 
  filter(id %in% c(baseline_magic, baseline_deepimpute)) %>%
  get_scores(syn, sub_df = .) %>%
  mutate(ensemble_label = all_sub_df$model_name[match(id, all_sub_df$id)]) %>%
  select(dataset, nrmse_score, spearman_score, ensemble_label) 

all_scores <- bind_rows(all_scores, baseline_scores)

all_scores$nrmse_score <- -all_scores$nrmse_score

# rank submissions all ensemble models and baseline
rank_df <- all_scores %>%
  rank_submissions("nrmse_score", "spearman_score", "ensemble_label") %>%
  mutate(team = gsub("\\d? \\+ (.*)$", "\\1", ensemble_label),
         ensemble_label = factor(ensemble_label, levels = c("0 + GOAL_LAB", 
                                                            "1 + DLS5", 
                                                            "2 + LDExplore",
                                                            "Baseline MAGIC",
                                                            "3 + Metformin-121",
                                                            "Baseline DeepImpute")))

# plot
# pal_simpsons(alpha = 0.9)(9)
my_labels <- c(primary_rank = "NRMSE", secondary_rank = "Spearman Correlation")
my_colors <- c(primary_rank = "#709AE1E5", secondary_rank = "#FED439E5")
rank_line_p <- ggplot(rank_df %>% gather("metrics", "ranks", c(2:3)),
       aes(ensemble_label, 1/ranks, color = metrics, group = metrics)) +
  geom_line(linewidth = 1) +
  # ggrepel::geom_label_repel(data = . %>% filter(overall_rank < 3, metrics == "primary_rank"),
  #                           aes(label = team), box.padding = unit(0.4, "lines"),
  #                           point.padding = unit(0.2, "lines")) +
  scale_color_manual(values = my_colors, labels = my_labels) +
  scale_x_discrete(guide = guide_axis(angle = 45)) + 
  labs(x = NULL, y = "1 / Ranks", color = "Metrics") +
  theme_bw(base_size = 18) + 
  theme(legend.position = c(0.95, 0.95), legend.justification = c("right", "top"),
        panel.grid.minor.y = element_blank(),
        text = element_text(size = 18),
        axis.title = element_text(size = 20))

  
## boot
# get corresponding metric names and for each task
metrics_lookup <- c("nrmse_score", "spearman_score")

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
  gather("metrics", "ranks", c(primary_rank, secondary_rank))

# plot
p_top1 <- bf_df %>%
  filter(metrics == "primary_rank") %>%
  mutate(groups = factor(
    case_when(
      primary_bf > 0 & primary_bf <= 5  ~ "< 5",
      primary_bf > 5 & primary_bf <= 30 ~ "5 - 30",
      primary_bf > 30 ~ "> 30",
      TRUE ~ "Reference"
    ), 
    levels = c("Reference", "< 5", "5 - 30", "> 30")),
  ensemble_label = factor(ensemble_label, levels = c("0 + GOAL_LAB", 
                                                     "1 + DLS5", 
                                                     "2 + LDExplore",
                                                     "Baseline MAGIC",
                                                     "3 + Metformin-121",
                                                     "Baseline DeepImpute"))) %>%
  ggplot(aes(ensemble_label, 1/ranks, color = groups)) + 
  labs(x = NULL, y = "1 / (Bootstrapped ranks of NRMSE)", color = "Bayes Factor") +
  geom_boxplot(lwd = 1.2, fatten = 1) + 
  scale_x_discrete(limits=rev) +
  coord_flip() +
  theme_classic(base_size = 16) + 
  scale_color_manual(values = c(
    "Reference" = "#A81A50", 
    '< 5' = '#F94551', 
    "5 - 30" = "#FCB335",
    "> 30" = "#32A0B5"
  ), drop = FALSE) +
  theme(text = element_text(size = 16),
        axis.title = element_text(size = 18))

p_top2 <- bf_df %>%
  filter(metrics == "secondary_rank") %>%
  mutate(groups = factor(
    case_when(
      secondary_bf > 0 & secondary_bf <= 5  ~ "< 5",
      secondary_bf > 5 & secondary_bf <= 30 ~ "5 - 30",
      secondary_bf > 30 ~ "> 30",
      TRUE ~ "Reference"
    ), 
    levels = c("Reference", "< 5", "5 - 30", "> 30")),
    ensemble_label = factor(ensemble_label, levels = c("0 + GOAL_LAB", 
                                                       "1 + DLS5", 
                                                       "2 + LDExplore",
                                                       "Baseline MAGIC",
                                                       "3 + Metformin-121",
                                                       "Baseline DeepImpute"))) %>%
  ggplot(aes(ensemble_label, 1/ranks, color = groups)) + 
  labs(x = NULL, y = "1 / (Bootstrapped Ranks of Spearman Correlation)", color = "Bayes Factor") +
  geom_boxplot(lwd = 1.2, fatten = 1) + 
  scale_x_discrete(limits=rev) +
  coord_flip() +
  theme_classic(base_size = 16) + 
  scale_color_manual(values = c(
    "Reference" = "#A81A50", 
    '< 5' = '#F94551', 
    "5 - 30" = "#FCB335",
    "> 30" = "#32A0B5"
  ), drop = FALSE) +
  theme(text = element_text(size = 16),
        axis.title = element_text(size = 18))

pdf(file="sc1_ensemble_analysis.pdf", width = 16, height = 10)
rank_line_p;
p_top1 / p_top2 + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "top", legend.direction = "horizontal")
dev.off()