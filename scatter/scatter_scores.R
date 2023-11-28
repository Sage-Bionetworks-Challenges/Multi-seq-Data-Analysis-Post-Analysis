source("utils/setup.R")
source("utils/bootstrap_funcs.R")
source("utils/synapse_funcs.R")
source("utils/plot_funcs.R")

task_n <- 1
metrics <- metrics_lookup[[task_n]]


# Reading submission data -------------------------------------------------
sub_data <- file.path(data_dir, str_glue("final_submissions_task{task_n}.rds"))
score_data <- file.path(data_dir, str_glue("final_scores_task{task_n}.rds"))
if (!all(file.exists(sub_data, score_data))) source("submission/get_submissions.R")

sub_df <- readRDS(sub_data)
scores_df <- readRDS(score_data)

# re-rank all submissions
sub_ranks <- rank_submissions(scores_df %>% mutate(!!sym(metrics[1]) := -!!sym(metrics[1])), metrics[1], metrics[2])
sub_df <- sub_df %>% arrange(match(team, sub_ranks$team))

# Label model names -------------------------------------------------------
baseline_magic <- "9732066"
baseline_deepimpute <- "9732074"
top_performer <- sub_df$id[1]
sub_df <- sub_df %>% mutate(model_name = case_when(
  id == baseline_magic ~ "Baseline MAGIC",
  id == baseline_deepimpute ~ "Baseline DeepImpute",
  TRUE ~ as.character(team)
))


# prepare df for plotting ----------------------------------------------------------
scores_df$model_name <- sub_df$model_name[match(scores_df$id, sub_df$id)]

sc1_plot_df <- scores_df %>%
  separate(dataset, into = c("prefix", "prop", "replicate"), sep = "_", remove = FALSE) %>%
  mutate(
    downsampled = if_else(prop %in% c("p00625", "p0125", "p025"), "by_cells", "by_reads"),
    prefix = toupper(prefix),
    prop = case_when(
      prop == "p00625" ~ "6.25%",
      prop == "p0125" ~ "12.5%",
      prop == "p025" ~ "25%",
      TRUE ~ gsub("p(\\d+k)", "\\1", prop)
    ),
    prop = forcats::fct_relevel(prop, "6.25%", "12.5%", "25%", "10k", "20k", "50k"),
    model_name = case_when(
      model_name == "zoradeng_post" ~ "Anonymous team 1",
      model_name == "moqri_post" ~ "Anonymous team 2",
      model_name == "USF biostat_post" ~ "Anonymous team 3",
      model_name == "BBKCS_post" ~ "BBKCS",
      TRUE ~ model_name
    ),
    model_name = forcats::fct_relevel(
      model_name, "GOAL_LAB", "DLS5",
      "Anonymous team 1", "BBKCS",
      "Anonymous team 2", "LDExplore",
      "Anonymous team 3", "Baseline MAGIC",
      "Baseline DeepImpute",
      "Metformin-121"
    )
  )



# plotting ----------------------------------------------------------
# downsampled by reads
sc1_p1 <- ggplot(
  sc1_plot_df %>% filter(downsampled == "by_reads"),
  aes(x = nrmse_score, y = spearman_score)
) +
  geom_jitter(aes(color = model_name), size = 2, alpha = 0.8) +
  ggsci::scale_color_jco() +
  theme_bw(base_size = 16) +
  facet_grid(prop ~ prefix, scales = "free") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme(strip.background = element_blank()) +
  labs(x = "NRMSE", y = "Spearman Correlation", color = "Model", title = "Metrics Scatter Plot - Task 1", subtitle = "Dowsampled by reads")

# downsampled by cells
sc1_p2 <- ggplot(
  sc1_plot_df %>% filter(downsampled == "by_cells"),
  aes(x = nrmse_score, y = spearman_score)
) +
  geom_jitter(aes(color = model_name), size = 2, alpha = 0.8) +
  ggsci::scale_color_jco() +
  theme_bw(base_size = 16) +
  facet_grid(prop ~ prefix) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme(strip.background = element_blank()) +
  labs(x = "NRMSE", y = NULL, color = "Model", subtitle = "Dowsampled by cells")

# combine
sc1_p <- sc1_p1 + sc1_p2 +
  plot_layout(guides = "collect", widths = c(5, 2)) &
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.box.background = element_rect(colour = "black"))

pdf(file = "metrics_scatter.pdf", width = 12, height = 6)
sc1_p
dev.off()


# Upload to synapse------------------------------------------------------------------
saveRDS(sc1_plot_df, "data/metrics_scatter_data.rds")
file <- synapseclient$File("data/metrics_scatter_data.rds", parent = "syn51270280")
stored <- syn$store(file, used = "syn51320982")
file <- synapseclient$File("metrics_scatter.pdf", parent = "syn51270280")
stored <- syn$store(file, used = "syn51320982")
