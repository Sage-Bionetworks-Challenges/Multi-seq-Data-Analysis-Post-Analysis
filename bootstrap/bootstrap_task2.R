source("utils/setup.R")

task_n <- 2


# Reading submission data -------------------------------------------------
sub_data <- file.path(data_dir, str_glue("final_submissions_task{task_n}.rds"))
score_data <- file.path(data_dir, str_glue("final_scores_task{task_n}.rds"))
if (!all(file.exists(sub_data, score_data))) source("submission/get_submissions.R")


sub_df <- readRDS(sub_data)
scores_df <- readRDS(score_data)

# label model names
baseline_macs2 <- "9732044"
top_performer <- sub_df$id[1]
sub_df <- sub_df %>% mutate(model_name = case_when(id == baseline_macs2 ~ "Baseline Macs2",
                                                   TRUE ~ as.character(team)))


# Retrieve scores ---------------------------------------------------------
scores_df <- get_scores(syn, sub_df)


# Bootstrapping -----------------------------------------------------------
# bootstrapping the rankings
metrics <- metrics_lookup[[task_n]]
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
# compute bayes factor between each model
ref_ids <- c(baseline_macs2, top_performer)
ref_names <- c("baseline_macs2", "top_performer")
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


# Plotting ----------------------------------------------------------------
# against top performer
p_top1 <- bf_df %>%
  filter(ref_model == "top_performer", metrics == "primary_rank") %>%
  bootstrap_boxplot(model_name, ranks, primary_bf, 
                    bf_cutoffs=c(5, 30), ref_label = "Ref: Top Performer") +
  labs(x = NULL, y = "1 / (Bootstrapped Ranks of Summed Scores)", color = "Bayes Factor")

p_top2 <- bf_df %>%
  filter(ref_model == "top_performer", metrics == "secondary_rank") %>%
  bootstrap_boxplot(model_name, ranks, primary_bf, 
                    bf_cutoffs=c(5, 30), ref_label = "Ref: Top Performer") +
  labs(x = NULL, y = "1 / (Bootstrapped Ranks of Jaccard similarity)", color = "Bayes Factor")

p_top <- p_top1 / p_top2 + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "top", legend.direction = "horizontal")

# against baseline MACS2
p_macs1 <- bf_df %>%
  filter(ref_model == "baseline_magic", metrics == "primary_rank") %>%
  bootstrap_boxplot(model_name, ranks, primary_bf, bf_cutoffs=c(5, 30)) +
  labs(x = NULL, y = "1 / (Bootstrapped Ranks of Summed Scores)", color = "Bayes Factor")

p_macs2 <- bf_df %>%
  filter(ref_model == "baseline_magic", metrics == "secondary_rank") %>%
  bootstrap_boxplot(model_name, ranks, primary_bf, bf_cutoffs=c(5, 30)) +
  labs(x = NULL, y = "1 / (Bootstrapped Ranks of Jaccard similarity)", color = "Bayes Factor")

p_macs <- p_macs1 / p_macs2 + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "top", legend.direction = "horizontal")

pdf(file="sc2_bootstrap_bayes_factor.pdf", width = 18, height = 12)
p_top; p_macs
dev.off()

