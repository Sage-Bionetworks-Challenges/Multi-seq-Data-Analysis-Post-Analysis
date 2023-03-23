library(ggplot2)
library(dplyr) # 1.1.0
library(tidyr) # 1.3.0
library(stringr)

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
library(patchwork)

# top performer
p_top1 <- bf_df %>%
  filter(ref_model == "top_performer", metrics == "primary_rank") %>%
  mutate(groups = as.factor(case_when(
    primary_bf < 5 & (id != top_performer)  ~ "Bayes Factor < 5",
    primary_bf >= 5 & (id != top_performer)  ~ "Bayes Factor >= 5",
    id == top_performer ~ "Top Performer"
  )),
  model_name = factor(model_name, levels = unique(sub_df$model_name))) %>% 
  ggplot(aes(model_name, 1/ranks, color = groups)) + 
  labs(title = str_glue("Bootstrapped submissions agsinst Top Performers"), 
       subtitle = "Summed Scores (jaccard_similarity + recall_ubiquitous + recall_tss [+ recall_cellspecific for mouse data])",
       x = NULL, y = NULL, color = NULL) +
  geom_boxplot(lwd = 1.2, fatten = 1) + 
  scale_x_discrete(limits=rev) +
  coord_flip() +
  theme_classic(base_size = 16) + 
  scale_color_manual(values = c(
    "Top Performer" = "#FE4365", 
    'Bayes Factor < 5' = '#FC9D9A', 
    "Bayes Factor >= 5" = "#C8C8A9"
  )) +
  theme(
    text = element_text(size = 16),
    plot.title = element_text(size=22),
    axis.text.x.bottom = element_text(size = 18),
    axis.title.y=element_text(size = 24), 
    axis.text.y=element_text(size = 18))

# 2rd metric
p_top2 <- bf_df %>%
  filter(ref_model == "top_performer", metrics == "secondary_rank") %>%
  mutate(groups = as.factor(case_when(
    secondary_bf < 5 & (id != top_performer)  ~ "Bayes Factor < 5",
    secondary_bf >= 5 & (id != top_performer) ~ "Bayes Factor >= 5",
    id == top_performer ~ "Top Performer"
    
  )),
  model_name = factor(model_name, levels = arrange(., metrics) %>% pull(model_name) %>% unique())) %>% 
  ggplot(aes(model_name, 1/ranks, color = groups)) + 
  labs(subtitle = "Jaccard similarity",
       x = NULL, y = "1 / Bootstrapped Ranks", color = NULL) +
  geom_boxplot(lwd = 1.2, fatten = 1) + 
  scale_x_discrete(limits=rev) +
  coord_flip() +
  theme_classic(base_size = 16) + 
  scale_color_manual(values = c(
    "Top Performer" = "#FE4365", 
    'Bayes Factor < 5' = '#FC9D9A', 
    "Bayes Factor >= 5" = "#C8C8A9")) +
  theme(
    text = element_text(size = 16),
    plot.title = element_text(size=22),
    axis.text.x.bottom = element_text(size = 18),
    axis.title.y=element_text(size = 24), 
    axis.text.y=element_text(size = 18))

p_top <- p_top1 / p_top2 + plot_layout(guides = "collect")

### Baseline macs2
p_macs2.1 <- bf_df %>%
  filter(ref_model == "baseline_macs2", metrics == "primary_rank") %>%
  mutate(groups = as.factor(case_when(
    primary_bf < 5 & (id != baseline_macs2)  ~ "Bayes Factor < 5",
    primary_bf >= 5 & (id != baseline_macs2)  ~ "Bayes Factor >= 5",
    id == baseline_macs2 ~ "Baseline Macs2"
  )),
  model_name = factor(model_name, levels = unique(sub_df$model_name))) %>% 
  ggplot(aes(model_name, 1/ranks, color = groups)) + 
  labs(title = str_glue("Bootstrapped submissions agsinst Baseline MAGIC"), 
       subtitle = "Summed Scores (jaccard_similarity + recall_ubiquitous + recall_tss [+ recall_cellspecific for mouse data])",
       x = NULL, y = NULL, color = NULL) +
  geom_boxplot(lwd = 1.2, fatten = 1) + 
  scale_x_discrete(limits=rev) +
  coord_flip() +
  theme_classic(base_size = 16) + 
  scale_color_manual(values = c(
    "Baseline Macs2" = "#FE4365", 
    'Bayes Factor < 5' = '#FC9D9A', 
    "Bayes Factor >= 5" = "#C8C8A9")) +
  theme(
    text = element_text(size = 16),
    plot.title = element_text(size=22),
    axis.text.x.bottom = element_text(size = 18),
    axis.title.y=element_text(size = 24), 
    axis.text.y=element_text(size = 18))

# 2rd metric
p_macs2.2 <- bf_df %>%
  filter(ref_model == "baseline_macs2", metrics == "secondary_rank") %>%
  mutate(groups = as.factor(case_when(
    secondary_bf < 5 & (id != baseline_macs2)  ~ "Bayes Factor < 5",
    secondary_bf >= 5 & (id != baseline_macs2) ~ "Bayes Factor >= 5",
    id == baseline_macs2 ~ "Baseline Macs2"
  )),
  model_name = factor(model_name, levels = arrange(., metrics) %>% pull(model_name) %>% unique())) %>% 
  ggplot(aes(model_name, 1/ranks, color = groups)) + 
  labs(subtitle = "Jaccard similarity",
       x = NULL, y = "1 / Bootstrapped Ranks", color = NULL) +
  geom_boxplot(lwd = 1.2, fatten = 1) + 
  scale_x_discrete(limits=rev) +
  coord_flip() +
  theme_classic(base_size = 16) + 
  scale_color_manual(values = c(
    "Baseline Macs2" = "#FE4365", 
    'Bayes Factor < 5' = '#FC9D9A', 
    "Bayes Factor >= 5" = "#C8C8A9")) +
  theme(
    text = element_text(size = 16),
    plot.title = element_text(size=22),
    axis.text.x.bottom = element_text(size = 18),
    axis.title.y=element_text(size = 24), 
    axis.text.y=element_text(size = 18))

p_macs2 <- p_macs2.1 / p_macs2.2 + plot_layout(guides = "collect")

pdf(file="sc2_bootstrap_bayes_factor.pdf", width = 18, height = 12)
p_top; p_macs2
dev.off()

