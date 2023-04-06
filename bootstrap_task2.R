library(ggplot2)
library(dplyr) # 1.1.0
library(tidyr) # 1.3.0
library(stringr)
library(patchwork)

source("utils.R")
ncores <- parallel::detectCores() - 1

# set up synapse
reticulate::use_condaenv('synapse')
synapseclient <- reticulate::import('synapseclient')
syn <- synapseclient$Synapse()
syn$login(silent = TRUE)

# download the submissions
view_id <- "syn51157023"
eval_id <- "9615024"
sub_df <- get_ranked_submissions(syn, view_id, eval_id, "private")

# label model names
baseline_macs2 <- "9732044"
top_performer <- sub_df$id[1]
sub_df <- sub_df %>% mutate(model_name = case_when(id == baseline_macs2 ~ "Baseline Macs2",
                                                   TRUE ~ as.character(team)))

# get all scores
scores_df <- get_scores(syn, sub_df)

# get corresponding metric names and for each task
metrics_lookup <- c("summed_score", "jaccard_similarity")

# bootstrapping the rankings
boot_df <- bootstrap(.data = scores_df,
                     seq_size = length(unique(scores_df$dataset)),
                     .by = "id",
                     n_iterations = 1000,
                     seed = 98109,
                     ncores = ncores)

# rank submissions for all bootstraps
rank_df <- boot_df %>%
  nest(scores = c(metrics_lookup[1], metrics_lookup[2], dataset, id, team), .by = bs_n) %>% 
  mutate(ranks = parallel::mclapply(scores, rank_submissions, metrics_lookup[1], metrics_lookup[2], mc.cores = ncores)) %>%
  unnest(cols = ranks) %>%
  select(-scores) %>%
  mutate(model_name = sub_df$model_name[match(id, sub_df$id)])

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

# plotting
# top performer
p_top1 <- bf_df %>%
  filter(ref_model == "top_performer", metrics == "primary_rank") %>%
  mutate(groups = factor(
    case_when(
      primary_bf > 0 & primary_bf <= 5  ~ "< 5",
      primary_bf > 5 & primary_bf <= 30 ~ "5 - 30",
      primary_bf > 30 ~ "> 30",
      TRUE ~ "Ref: Top Performer"
    ), 
    levels = c("Ref: Top Performer", "< 5", "5 - 30", "> 30")),
  model_name = factor(model_name, levels = unique(sub_df$model_name))) %>% 
  ggplot(aes(model_name, 1/ranks, color = groups)) + 
  labs(x = NULL, y = "1 / (Bootstrapped Ranks of Summed Scores)", color = "Bayes Factor") +
  geom_boxplot(lwd = 1.2, fatten = 1) + 
  scale_x_discrete(limits=rev) +
  coord_flip() +
  theme_classic(base_size = 16) + 
  scale_color_manual(values = c(
    "Ref: Top Performer" = "#A81A50", 
    '< 5' = '#F94551', 
    "5 - 30" = "#FCB335",
    "> 30" = "#32A0B5"
  ), drop = FALSE) +
  theme(text = element_text(size = 16),
        axis.title = element_text(size = 18))

# 2rd metric
p_top2 <- bf_df %>%
  filter(ref_model == "top_performer", metrics == "secondary_rank") %>%
  mutate(groups = factor(
    case_when(
      secondary_bf > 0 & secondary_bf <= 5  ~ "< 5",
      secondary_bf > 5 & secondary_bf <= 30 ~ "5 - 30",
      secondary_bf > 30 ~ "> 30",
      TRUE ~ "Ref: Top Performer"
    ), 
    levels = c("Ref: Top Performer", "< 5", "5 - 30", "> 30")),
  model_name = factor(model_name, levels = arrange(., metrics) %>% pull(model_name) %>% unique())) %>% 
  ggplot(aes(model_name, 1/ranks, color = groups)) + 
  labs(x = NULL, y = "1 / (Bootstrapped Ranks of Jaccard similarity)", color = "Bayes Factor") +
    geom_boxplot(lwd = 1.2, fatten = 1) + 
    scale_x_discrete(limits=rev) +
    coord_flip() +
    theme_classic(base_size = 16) + 
    scale_color_manual(values = c(
      "Ref: Top Performer" = "#A81A50", 
      '< 5' = '#F94551', 
      "5 - 30" = "#FCB335",
      "> 30" = "#32A0B5"
    ), drop = FALSE) +
    theme(text = element_text(size = 16),
          axis.title = element_text(size = 18))

p_top <- p_top1 / p_top2 + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "top", legend.direction = "horizontal")

### Baseline macs2
p_macs2.1 <- bf_df %>%
  filter(ref_model == "baseline_macs2", metrics == "primary_rank") %>%
  mutate(groups = factor(
    case_when(
      primary_bf > 0 & primary_bf <= 5  ~ "< 5",
      primary_bf > 5 & primary_bf <= 30 ~ "5 - 30",
      primary_bf > 30 ~ "> 30",
      TRUE ~ "Reference"
    ), 
    levels = c("Reference", "< 5", "5 - 30", "> 30")),
  model_name = factor(model_name, levels = unique(sub_df$model_name))) %>% 
  ggplot(aes(model_name, 1/ranks, color = groups)) + 
  labs(x = NULL, y = "1 / (Bootstrapped Ranks of Summed Scores)", color = "Bayes Factor") +
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

# 2rd metric
p_macs2.2 <- bf_df %>%
  filter(ref_model == "baseline_macs2", metrics == "secondary_rank") %>%
  mutate(groups = factor(
    case_when(
      secondary_bf > 0 & secondary_bf <= 5  ~ "< 5",
      secondary_bf > 5 & secondary_bf <= 30 ~ "5 - 30",
      secondary_bf > 30 ~ "> 30",
      TRUE ~ "Reference"
    ), 
    levels = c("Reference", "< 5", "5 - 30", "> 30")),
  model_name = factor(model_name, levels = arrange(., metrics) %>% pull(model_name) %>% unique())) %>% 
  ggplot(aes(model_name, 1/ranks, color = groups)) + 
  labs(x = NULL, y = "1 / (Bootstrapped Ranks of Jaccard similarity)", color = "Bayes Factor") +
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

p_macs2 <- p_macs2.1 / p_macs2.2 + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "top", legend.direction = "horizontal")

pdf(file="sc2_bootstrap_bayes_factor.pdf", width = 18, height = 12)
p_top; p_macs2
dev.off()

