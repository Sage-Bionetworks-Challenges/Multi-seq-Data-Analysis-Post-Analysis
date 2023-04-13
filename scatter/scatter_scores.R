library(ggplot2)
library(dplyr) 
library(tidyr)
library(stringr)
library(ggsci)
library(patchwork)

## Task 1 - Final Round Only ------------------------------------------------

# set up syanpse ----------------------------------------------------------
reticulate::use_condaenv('synapse')
synapseclient <- reticulate::import('synapseclient')
syn <- synapseclient$Synapse()
syn$login(silent = TRUE)


# download the submissions ------------------------------------------------
view_id <- "syn51157023"
eval_id <- "9615023"
sub_df <- get_ranked_submissions(syn, view_id, eval_id, "private")


# label model names -------------------------------------------------------
baseline_magic <- "9732066"
baseline_deepimpute <- "9732074"
top_performer <- sub_df$id[1]
sub_df <- sub_df %>% mutate(model_name = case_when(id == baseline_magic ~ "Baseline MAGIC",
                                                   id == baseline_deepimpute ~ "Baseline DeepImpute",
                                                   TRUE ~ as.character(team)))


# get all scores ----------------------------------------------------------
scores_df <- get_scores(syn, sub_df)
scores_df$model_name <- sub_df$model_name[match(scores_df$id, sub_df$id)]

sc1_plot_df <- scores_df %>%
  separate(dataset, into = c("prefix", "prop", "replicate"), sep = "_", remove = FALSE) %>%
  mutate(downsampled = if_else(prop %in% c("p00625", "p0125", "p025"), "by_cells", "by_reads"),
         model_name = factor(model_name, levels = unique(model_name)),
         prefix = toupper(prefix),
         prop = factor(case_when(prop == "p00625" ~ "6.25%",
                          prop == "p0125" ~ "12.5%",
                          prop == "p025" ~ "25%",
                          TRUE ~ gsub("p(\\d+k)", "\\1", prop)),
                       levels = c("6.25%", "12.5%", "25%", "10k", "20k", "50k")))


# plotting ----------------------------------------------------------
# downsampled by reads
sc1_p1 <- ggplot(sc1_plot_df %>% filter(downsampled == "by_reads"),
                 aes(x=nrmse_score, y=spearman_score)) +
  geom_jitter(aes(color = model_name), size = 2, alpha = 0.8) +
  scale_color_jco() +
  theme_bw(base_size = 16) +
  facet_grid(prop~prefix, scales = "free") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme(strip.background = element_blank()) +
  labs(x = "NRMSE", y = "Spearman Correlation", color = "Model", title = "Metrics Scatter Plot - Task 1", subtitle = "Dowsampled by reads")

# downsampled by cells
sc1_p2 <- ggplot(sc1_plot_df %>% filter(downsampled == "by_cells"),
                 aes(x=nrmse_score, y=spearman_score)) +
  geom_jitter(aes(color = model_name), size = 2, alpha = 0.8) +
  scale_color_jco() +
  theme_bw(base_size = 16) +
  facet_grid(prop~prefix) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme(strip.background = element_blank()) +
  labs(x = "NRMSE", y = NULL, color = "Model", subtitle = "Dowsampled by cells") 

# combine
sc1_p <- sc1_p1 + sc1_p2 +
  plot_layout(guides = "collect", widths = c(5, 2))  &
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.box.background = element_rect(colour = "black"))


# visualize scores on each gene ------------------------------------------------
# read scores - generated via get_gene_scores.R
if (F) {
  gene_scores <- readRDS("all_genes_scores_sc1.rds")
  
  # plotting
  # downsampled by reads
  sc1_p1 <- ggplot(gene_scores %>% filter(downsampled == "by_reads"),
                   aes(x=nrmse_score, y=spearman_score)) +
    geom_point(aes(color = model_name), size = 0.8, alpha = 0.5) +
    scale_color_jco() +
    theme_bw(base_size = 16) +
    facet_grid(prop~prefix, scales = "free") +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    theme(strip.background = element_blank()) +
    labs(x = "NRMSE", y = "Spearman Correlation", color = "Model", title = "Metrics Scatter Plot All Genes - Task 1", subtitle = "Dowsampled by reads")
  
  # downsampled by cells
  sc1_p2 <- ggplot(gene_scores %>% filter(downsampled == "by_cells"),
                   aes(x=nrmse_score, y=spearman_score)) +
    geom_point(aes(color = model_name), size = 0.8, alpha = 0.5, hape='.') +
    scale_color_jco() +
    theme_bw(base_size = 16) +
    facet_grid(prop~prefix) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    theme(strip.background = element_blank()) +
    labs(x = "NRMSE", y = NULL, color = "Model", subtitle = "Dowsampled by cells") 
  
  # combine
  sc1_genes_p <- sc1_p1 + sc1_p2 +
    plot_layout(guides = "collect", widths = c(5, 2))  &
    theme(legend.position = "bottom", legend.direction = "horizontal", legend.box.background = element_rect(colour = "black"))
  
  
  pdf("all_gene_scores_sc1.pdf", width = 20, height = 12)
  sc1_genes_p
  dev.off()
}

## Task 2 - Final Round Only ------------------------------------------------
# repeat
view_id <- "syn51157023"
eval_id <- "9615024"
sub_df <- get_ranked_submissions(syn, view_id, eval_id, "private")

baseline_macs2 <- "9732044"
top_performer <- sub_df$id[1]
sub_df <- sub_df %>% mutate(model_name = case_when(id == baseline_macs2 ~ "Baseline MACS2",
                                                   TRUE ~ as.character(team)))
scores_df <- get_scores(syn, sub_df)
scores_df$model_name <- sub_df$model_name[match(scores_df$id, sub_df$id)]

sc2_p1 <- ggplot(scores_df %>%
       filter(grepl("ds1", dataset)) %>%
       separate(dataset, into =c("prefix", "dummy", "cluster", "replicate", "prop"), sep = "\\.", remove = FALSE) %>%
       mutate(model_name = factor(model_name, levels = unique(model_name)),
              prefix = toupper(prefix),
              prop = factor(paste0(as.numeric(sub(".*pct_(\\d+).bed", "\\1", dataset)), "%"),
                            levels = c( "1%", "5%", "10%"))),
       aes(x = summed_score, y = jaccard_similarity, color = model_name)) +
  geom_jitter(aes(color = model_name), size = 2, alpha = 0.8) +
  scale_color_jco() +
  theme_bw(base_size = 16) +
  facet_grid(prop~prefix) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme(strip.background = element_blank()) +
  labs(x = "Summed Score", y = "Jaccard Similarity", color = "Model", title = "Metrics Scatter Plot - Task 2")


sc2_p2 <- ggplot(scores_df %>%
       filter(grepl("ds2", dataset)) %>%
       separate(dataset, into =c("prefix", "dummy", "cluster", "replicate", "prop"), sep = "\\.", remove = FALSE) %>%
       mutate(model_name = factor(model_name, levels = unique(model_name)),
              prefix = toupper(prefix),
              prop = factor(paste0(as.numeric(sub(".*pct_(\\d+).bed", "\\1", dataset)), "%"),
                            levels = c( "10%", "20%", "50%"))),
       aes(x=summed_score, y=jaccard_similarity, color = model_name)) +
  geom_jitter(aes(color = model_name), size = 2, alpha = 0.8) +
  scale_color_jco() +
  theme_bw(base_size = 16) +
  facet_grid(prop~prefix) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme(strip.background = element_blank()) +
  labs(x = "Summed Score", y = "Jaccard Similarity", color = "Model")

sc2_p <- sc2_p1 + sc2_p2 + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.box.background = element_rect(colour = "black"))

pdf(file="metrics_scatter.pdf", width = 12, height = 6)
sc1_p; sc2_p
dev.off()


# upload to synapse------------------------------------------------------------------

file <- synapseclient$File("metrics_scatter.pdf", parent="syn51270280")
stored <-syn$store(file, used = "syn51320982")

# large, use screenshot instead
# file <- synapseclient$File("metrics_scatter_all_gens_sc1.pdf", parent="syn51270280")
# stored <-syn$store(file, used = "syn51320982")


# ## Task 1 - Both Rounds ------------------------------------------------
# 
# # download the submissions from post-submission table
# view_id <- "syn51157023"
# eval_id <- "9615324"
# query <- stringr::str_glue(
#   "
#   SELECT
#     re_submitterid, 
#   	re_id, 
#   	prediction_fileid,
#     submission_scores,
#   	re_submission_phase,
#   	re_overall_rank
#     FROM syn51157023
#     WHERE submission_status = 'SCORED'
#       AND status = 'ACCEPTED'
#       AND evaluationid = '9615324' OR evaluationid = '9615023'
#       AND re_id is NOT NULL
#     ORDER BY re_overall_rank
#     "
# )
# 
# # download the submissions ordered by overall rank
# sub_df <- syn$tableQuery(query)$asDataFrame() %>%
#   mutate(across(everything(), as.character),
#          team = as.character(sapply(re_submitterid, get_name, syn = syn))) %>%
#   setnames(c("submitterid", "id", "prediction_fileid", "submission_scores", "submission_phase", "ranks", "team")) %>%
#   mutate(id = ifelse(submission_phase == "public", submitterid, id)) %>%
#   tibble::remove_rownames()
# 
# # label model name - repeat steps from task1 - final round only
# baseline_magic <- "9732066"
# baseline_deepimpute <- "9732074"
# top_performer <- sub_df$id[1]
# sub_df <- sub_df %>% mutate(model_name = case_when(id == baseline_magic ~ "Baseline MAGIC",
#                                                    id == baseline_deepimpute ~ "Baseline DeepImpute",
#                                                    TRUE ~ as.character(team)))
# 
# 
# # get all scores ----------------------------------------------------------
# scores_df <- get_scores(syn, sub_df) %>%
#   mutate(model_name = sub_df$model_name[match(id, sub_df$id)],
#          phase = sub_df$submission_phase[match(id, sub_df$id)]) %>%
#   mutate(tag = paste0(model_name, "_", phase))
# 
# sc1_plot_df <- scores_df %>%
#   separate(dataset, into = c("prefix", "prop", "replicate"), sep = "_", remove = FALSE) %>%
#   mutate(downsampled = if_else(prop %in% c("p00625", "p0125", "p025"), "by_cells", "by_reads"),
#          tag = factor(tag, levels = unique(scores_df$tag)),
#          prefix = toupper(prefix),
#          prop = factor(case_when(prop == "p00625" ~ "6.25%",
#                                  prop == "p0125" ~ "12.5%",
#                                  prop == "p025" ~ "25%",
#                                  TRUE ~ gsub("p(\\d+k)", "\\1", prop)),
#                        levels = c("6.25%", "12.5%", "25%", "10k", "20k", "50k")))
# 
# 
# # plotting ----------------------------------------------------------
# # downsampled by reads
# cols <- c(ggsci::pal_jco()(7), ggsci::pal_flatui()(10)[3:8])
# sc1_p1 <- ggplot(sc1_plot_df %>% filter(downsampled == "by_reads"),
#                  aes(x=-nrmse_score, y=spearman_score)) +
#   geom_jitter(aes(color = tag), size = 2, alpha = 1) +
#   scale_color_manual(values=cols) +
#   theme_bw(base_size = 16) +
#   facet_grid(prop~prefix, scales = "free") +
#   theme(strip.background = element_blank()) +
#   labs(x = "NRMSE", y = "Spearman Correlation", color = "Model", title = "Metrics Scatter Plot - Task 1", subtitle = "Dowsampled by reads")
# 
# # downsampled by cells
# sc1_p2 <- ggplot(sc1_plot_df %>% filter(downsampled == "by_cells"),
#                  aes(x=nrmse_score, y=spearman_score)) +
#   geom_jitter(aes(color = model_name), size = 2, alpha = 1) +
#   scale_color_jco() +
#   theme_bw(base_size = 16) +
#   facet_grid(prop~prefix) +
#   theme(strip.background = element_blank()) +
#   labs(x = "NRMSE", y = NULL, color = "Model", subtitle = "Dowsampled by cells") 
# 
# # combine
# sc1_p_all <- sc1_p1 + sc1_p2 +
#   plot_layout(guides = "collect", widths = c(5, 2))  &
#   theme(legend.position = "bottom", legend.direction = "horizontal", legend.box.background = element_rect(colour = "black"))
# 
# 
# pdf(file="metrics_scatter_both_rounds_task1.pdf", width = 12, height = 6)
# sc1_p_all
# dev.off()
# 
# 
# # upload to synapse------------------------------------------------------------------
# 
# file <- synapseclient$File("metrics_scatter_both_rounds_task1.pdf", parent="syn51270280")
# # stored <-syn$store(file, used = "")

