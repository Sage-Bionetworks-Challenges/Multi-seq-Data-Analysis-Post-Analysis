library(ggplot2)
library(dplyr)
library(tidyr)
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
task_eval_id <- "9615023"
sub_df <- get_submissions(syn, view_id, task_eval_id)

# label model names
baseline_magic <- "9732066"
baseline_deepimpute <- "9732074"
top_performer <- sub_df$id[1]
sub_df <- sub_df %>%
  mutate(model_name = case_when(
    id == top_performer ~ "Top Performer",
    id == baseline_magic ~ "Baseline_MAGIC",
    id == baseline_deepimpute ~ "Baseline_DeepImpute",
    TRUE ~ str_glue("{team}_{id}")
  )
)

# get all scores
scores_df <- get_scores(syn, sub_df)

# get corresponding metric names and for each task
metrics_lookup <- c("nrmse_score", "spearman_score")
#   metrics_lookup <- c("summed_score", "jaccard_similarity")

# correct the direction of nrmse
scores_df$nrmse_score <- -scores_df$nrmse_score

# bootstrapping the rankings
boot_df <- bootstrap(.data = scores_df,
                      seq_size = length(unique(scores_df$dataset)),
                      .by = "id",
                      n_iterations = 2000,
                      seed = 98109,
                      ncores = ncores)

# rank submissions for all bootstraps
rank_df <- boot_df %>%
  nest(scores = c(metrics_lookup[1], metrics_lookup[2], dataset, id, team), .by = n_bs) %>%
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
    mutate(primary_bf = bayes_factor(primary_rank[id != ref_ids[i]], ref_ranks$primary_rank),
           secondary_bf = bayes_factor(secondary_rank[id != ref_ids[i]], ref_ranks$secondary_rank),
           ref_model = ref_names[i])
}) %>% 
  bind_rows() %>%
  gather("metrics", "ranks", c(primary_rank, secondary_rank)) %>%
  mutate(model_name = factor(model_name, levels = sub_df$model_name))

# more plots
bf_df %>%
  filter(ref_model == "top_performer", metrics == "primary_rank") %>%
  mutate(groups = as.factor(case_when(
         ranks < 5 && !id %in% ref_ids  ~ "Bayes Factor < 5",
         ranks >= 5 && !id %in% ref_ids ~ "Bayes Factor >= 5",
         id == top_performer ~ "Top Performer",
         id == baseline_magic ~ "Baseline_MAGIC",
         id == baseline_deepimpute ~ "Baseline_DeepImpute"
        ))) %>% 
  ggplot(aes(model_name, 1/ranks, color = groups)) + 
  labs(y = "1 / Bootstrapped Ranks", color = NULL) +
  geom_boxplot(lwd = 1.2, fatten = 1) + 
  scale_x_discrete(limits=rev) +
  coord_flip() +
  theme_classic() + 
  scale_color_manual(values = c(
    "Top Performer" = "#12F10E", 
    'Bayes Factor < 5' = '#2F88FB', 
    "Bayes Factor >= 5" = "#BFEAE6",
    "Baseline_MAGIC" = "#FC7E5A",
    "Baseline_DeepImpute" = "#FC7E5A")) +
  theme(
    text = element_text(size = 16),
    axis.text.x.bottom = element_text(size = 18),
    axis.title.y=element_text(size = 24), 
    axis.text.y=element_text(size = 18))

ggsave(
  file="sc1-bf.pdf",
  # plot=gridExtra::grid.arrange(plot.auc, plot.bayes.top, plot.bayes.ref, ncol = 3, widths = c(3, 1, 1)),
  width = 24,
  height = 12
)


