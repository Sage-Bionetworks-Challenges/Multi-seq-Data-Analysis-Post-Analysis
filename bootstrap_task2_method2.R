# download the submissions
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(bedr)

source("metrics.R")
source("utils.R")
ncores <- parallel::detectCores() - 1

# set up synapse
reticulate::use_condaenv('synapse')
synapseclient <- reticulate::import('synapseclient')
syn <- synapseclient$Synapse()
syn$login(silent = TRUE)


# download the submissions
task_n <- "task2"
view_id <- "syn51157023"
eval_id <- "9615024"
gs_id <- "syn35294386"
sub_phase <- "private"

# gs_ids <- list(task1 = "syn34612394", task2 = "syn35294386")
# each team should be ranked only on the best submission, aka one submission per team from this table
sub_df <- get_ranked_submissions(syn, view_id, eval_id, "private")

# download the ground truth
gs_path <- syn$get(gs_id)["path"]
all_gs <- readRDS(gs_path)

# iterate each submission to perform bootstrapping
all_scores <- tibble()

for (sub_n in 1:nrow(sub_df)) { # use for loop to save computation for bootstrapping
  
  sub_id <- sub_df$id[sub_n]
  pred_id <- sub_df$prediction_fileid[sub_n]
  team <- sub_df$team[sub_n]

  message(str_glue("Scoring {team} ..."))
  # download the submission's prediction file
  # it will take some time, since each prediction file is a large tarball
  pred_path <- syn$get(pred_id)$path
  
  # create a temp dir to store all prediction files for one submission
  temp_dir <- str_glue("{team}_{sub_id}")
  dir.create(temp_dir, showWarnings = FALSE)
  
  # decompress the prediction files into the temp dir
  message("Untaring submission to ", temp_dir, " ...")
  untar(pred_path, exdir = temp_dir)
  
  pred_files <- list.files(file.path(temp_dir, "output"), pattern = ".bed")
  
  # collect all bootstrapping scores for each submission
  scores_df <- lapply(seq_along(pred_files), function(i) {
    
    pred_file <- pred_files[i]
    # detect file prefix used to read gs
    info <- strsplit(pred_file, "\\.")[[1]]
    prefix <- info[1]
    jaccards <- list()

    # read gs
    gs <- all_gs$gs_data[[sub_phase]][[prefix]]
    gs_ranked_filtered <- gs[["gs_ranked_filtered"]]
    gs_sort <- gs[["gs_sort"]]
    gs_sort_ubi <- gs[["gs_sort_ubi"]]
    gs_sort_tss <- gs[["gs_sort_tss"]]
    gs_sort_csp <- gs[["gs_sort_csp"]]
  
    # read prediction
    pred <- data.table::fread(file.path(temp_dir, "output", pred_file), data.table = FALSE, verbose = FALSE)

    jaccard_similarity <- 0
    recall_ubiquitous <- 0
    recall_tss <- 0
    recall_cellspecific <- 0
    summed_score <- 0

    # bootstrap 1000 times
    # make indices of each file are the same for all submissions
    bs_n <- 1000
    bs_indices <- boot_indices(seq_size = nrow(pred), n_iterations = bs_n, seed = 98109)

    bs_scores <- pbmcapply::pbmclapply(seq_along(bs_indices), function(n) {
      
      if (nrow(pred) > 0) {
        bs_indice <- bs_indices[[n]]

        colnames(pred)[1:3] <- c("chr", "start", "end")
        pred <- pred %>% filter(grepl("chr(\\d|X|Y)", chr))
        pred_sort <- bedr.merge.region(bedr.sort.region(data.frame(pred)[bs_indice, 1:3], check.zero.based = FALSE, check.chr = FALSE, check.valid = FALSE, verbose = FALSE), verbose = FALSE)

        # caculate scores
        j_result <- jaccard(gs_sort, pred_sort, check.merge = TRUE, check.chr = FALSE, check.sort = FALSE, check.valid = FALSE, verbose = FALSE)
        jaccard_similarity <- as.numeric(j_result$jaccard)
        recall_ubiquitous <- category_recall(gs_sort_ubi, pred_sort)
        recall_tss <- category_recall(gs_sort_tss, pred_sort)
        recall_cellspecific <- category_recall(gs_sort_csp, pred_sort)

        summed_score <- jaccard_similarity + recall_ubiquitous + recall_tss

        # report cellspecific for ds2, but not added to summed score
        if (prefix == "ds1") summed_score <- summed_score + recall_cellspecific
      }

      print(str_glue("Done {i} / {length(pred_files)}"))
      # collect scores for each test case
      return(
        tibble(
          dataset = pred_file,
          jaccard_similarity = jaccard_similarity,
          recall_ubiquitous = recall_ubiquitous,
          recall_tss = recall_tss,
          recall_cellspecific = recall_cellspecific,
          summed_score = summed_score,
          bs_n = n
        )
      )
    }, mc.cores = ncores) %>% bind_rows()
  }) %>% 
    bind_rows() %>%
    mutate(id = sub_id, team = team)
  
  # remove scored submission prediction files and temp dir to clean space
  unlink(temp_dir, recursive = TRUE)
  # append the each submission's bootstrapping results
  all_scores <- rbind(all_scores, scores_df)
}

# get corresponding metric names and for each task
metrics_lookup <- c("summed_score", "jaccard_similarity")

# rank submissions for all bootstraps
rank_df <- all_scores %>%
  nest(scores = c(metrics_lookup[1], metrics_lookup[2], dataset, id, team), .by = bs_n) %>%
  mutate(ranks = parallel::mclapply(scores, rank_submissions, metrics_lookup[1], metrics_lookup[2], mc.cores = ncores)) %>%
  unnest(cols = ranks) %>%
  select(-scores) %>%
  mutate(model_name = sub_df$model_name[match(id, sub_df$id)])

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
  file="sc2-bf.pdf",
  # plot=gridExtra::grid.arrange(plot.auc, plot.bayes.top, plot.bayes.ref, ncol = 3, widths = c(3, 1, 1)),
  width = 24,
  height = 12
)