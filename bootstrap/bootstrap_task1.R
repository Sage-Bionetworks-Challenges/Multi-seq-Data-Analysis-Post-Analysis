source("utils/setup.R")
source("utils/synapse_funcs.R")
source("utils/bootstrap_funcs.R")
source("utils/plot_funcs.R")

task_n <- 1
metrics <- metrics_lookup[[task_n]]


# Reading submission data -------------------------------------------------
sub_data <- file.path(data_dir, str_glue("final_submissions_task{task_n}.rds"))
score_data <- file.path(data_dir, str_glue("final_scores_task{task_n}.rds"))
if (!all(file.exists(sub_data, score_data))) source("submission/get_submissions.R")

sub_df <- readRDS(sub_data)
scores_df <- readRDS(score_data)
# correct the direction of nrmse
scores_df$nrmse_score <- -scores_df$nrmse_score

# re-rank all data
sub_ranks <- rank_submissions(scores_df, metrics[1], metrics[2])
sub_df <- sub_df %>% arrange(match(team, sub_ranks$team))

# label model names
baseline_magic <- "9732066"
baseline_deepimpute <- "9732074"
top_performer <- sub_df$id[1]
sub_df <- sub_df %>% 
  mutate(model_name = case_when(id == baseline_magic ~ "Baseline MAGIC",
                                id == baseline_deepimpute ~ "Baseline DeepImpute",
                                team == "zoradeng_post" ~ "Anonymous team 1",
                                team == "moqri_post" ~ "Anonymous team 2",
                                team == "USF biostat_post" ~ "Anonymous team 3",
                                TRUE ~ as.character(team)))


# Bootstrapping -----------------------------------------------------------
# bootstrapping the test case's scores
boot_df <- simple_bootstrap(.data = scores_df,
                            seq_size = length(unique(scores_df$dataset)),
                            .by = "id",
                            n_iter = 1000,
                            seed = 98109,
                            ncores = ncores)

# rank submissions for all bootstraps
rank_df <- boot_df %>%
  nest(scores = c(metrics[1], metrics[2], dataset, id, team), .by = bs_n) %>%
  mutate(ranks = parallel::mclapply(scores, rank_submissions, metrics[1], metrics[2], mc.cores = ncores)) %>%
  unnest(cols = ranks) %>%
  select(-scores) %>%
  mutate(model_name = sub_df$model_name[match(id, sub_df$id)])


# Compute bayes factors ---------------------------------------------------
# compute bayes factor of each model against reference
ref_ids <- c(baseline_magic, baseline_deepimpute, top_performer)
ref_names <- c("baseline_magic", "baseline_deepimpute", "top_performer")
bf_df <- lapply(seq_along(ref_ids), function(i) {
  ref_ranks <- rank_df %>% filter(id == ref_ids[i])
  rank_df %>%
    group_by(id) %>%
    mutate(primary_bf = calculate_bf(primary_rank, ref_ranks$primary_rank),
           secondary_bf = calculate_bf(secondary_rank, ref_ranks$secondary_rank),
           ref_model = ref_names[i])
  }) %>% 
  bind_rows() %>%
  gather("metrics", "ranks", c(primary_rank, secondary_rank)) %>%
  mutate(model_name = factor(model_name, levels = unique(sub_df$model_name)))

# Plotting ----------------------------------------------------------------
# against top performe
p_top1 <- bf_df %>%
  filter(ref_model == "top_performer", metrics == "primary_rank") %>%
  bootstrap_boxplot(model_name, ranks, primary_bf, ref = "GOAL_LAB",
                    bf_cutoffs = 3, ref_label = "Ref") +
  labs(x = NULL, y = "Bootstrapped ranks of NRMSE", color = NULL) +
  theme(legend.position = c(0.95, 0.85), legend.justification = c("right", "top"),
        legend.background = element_rect(color = "black", fill = "white", linewidth = 0.1),
        panel.grid.minor.y = element_blank()) 
# +
#   geom_vline(xintercept = 9.5, color = "#F94551", expand = c(0.1, 0.1), linetype = "dashed")

p_top2 <- bf_df %>%
  filter(ref_model == "top_performer", metrics == "secondary_rank") %>%
  bootstrap_boxplot(model_name, ranks, primary_bf, ref = "GOAL_LAB",
                    bf_cutoffs = 3, ref_label = "Ref") +
  labs(x = NULL, y = "Bootstrapped ranks of Spearman correlation", color = NULL) +
  theme(legend.position = c(0.95, 0.85), legend.justification = c("right", "top"),
        legend.background = element_rect(color = "black", fill = "white", linewidth = 0.1),
        panel.grid.minor.y = element_blank())
  
p_top <- p_top1 / p_top2

# against baseline MAGIC
p_magic1 <- bf_df %>%
  filter(ref_model == "baseline_magic", metrics == "primary_rank") %>%
  bootstrap_boxplot(model_name, ranks, primary_bf, ref = "Baseline MAGIC", ref_label = "Ref", bf_cutoffs = 3) +
  labs(x = NULL, y = "Bootstrapped ranks of NRMSE", color = NULL) +
  theme(legend.position = c(0.95, 0.85), legend.justification = c("right", "top"),
        legend.background = element_rect(color = "black", fill = "white", linewidth = 0.1),
        panel.grid.minor.y = element_blank())
  
p_magic2 <- bf_df %>%
  filter(ref_model == "baseline_magic", metrics == "secondary_rank") %>%
  bootstrap_boxplot(model_name, ranks, primary_bf,ref = "Baseline MAGIC", ref_label = "Ref", bf_cutoffs = 3) +
  labs(x = NULL, y = "Bootstrapped ranks of Spearman correlation", color = NULL) +
  theme(legend.position = c(0.95, 0.85), legend.justification = c("right", "top"),
        legend.background = element_rect(color = "black", fill = "white", linewidth = 0.1),
        panel.grid.minor.y = element_blank())

p_magic <- p_magic1 / p_magic2

# against baseline DeepImpute
p_dp1 <- bf_df %>%
  filter(ref_model == "baseline_deepimpute", metrics == "primary_rank") %>%
  bootstrap_boxplot(model_name, ranks, primary_bf, ref = "Baseline DeepImpute", ref_label = "Ref", bf_cutoffs = 3) +
  labs(x = NULL, y = "Bootstrapped ranks of NRMSE", color = NULL) +
  theme(legend.position = c(0.95, 0.85), legend.justification = c("right", "top"),
        legend.background = element_rect(color = "black", fill = "white", linewidth = 0.1),
        panel.grid.minor.y = element_blank())
  
p_dp2 <- bf_df %>%
  filter(ref_model == "baseline_deepimpute", metrics == "secondary_rank") %>%
  bootstrap_boxplot(model_name, ranks, primary_bf, ref = "Baseline DeepImpute", ref_label = "Ref", bf_cutoffs = 3) +
  labs(x = NULL, y = "Bootstrapped ranks of Spearman correlation", color = NULL) +
  theme(legend.position = c(0.95, 0.85), legend.justification = c("right", "top"),
        legend.background = element_rect(color = "black", fill = "white", linewidth = 0.1),
        panel.grid.minor.y = element_blank())
  
p_dp <- p_dp1 / p_dp2
  
# pdf(file="sc1_bootstrap_bayes_factor.pdf", width = 18, height = 12)
pdf(file="sc1_bootstrap_bayes_factor_post.pdf", width = 18, height = 12)
p_top; p_magic; p_dp
dev.off()

