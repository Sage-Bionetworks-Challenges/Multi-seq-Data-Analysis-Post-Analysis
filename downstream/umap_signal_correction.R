library(ggplot2)
library(dplyr)
library(stringr)

umap_data <- readRDS(file.path("data", "downstream_results_summary.rds"))$umap
ds3_umap_data <- umap_data %>% 
  filter(dataset == "ds3", downsampled_by == "reads") %>%
  mutate(
    proportion =  forcats::fct_relevel(proportion, "10k", "20k", "50k"),
    team = case_when(team == "zoradeng" ~ "Anon. team 1",
                     team == "moqri" ~ "Anon. team 2",
                     team == "USF-biostat" ~ "Anon. team 3",
                     team == "Baseline_MAGIC" ~ "MAGIC",
                     team == "Baseline_DeepImpute" ~ "DeepImpute",
                     TRUE ~ team),
    team = forcats::fct_relevel(team, "no_imp", "GOAL_LAB", "DLS5", "Anon. team 1",
                                "BBKCS", "Anon. team 2", "LDExplore",
                                "Anon. team 3", "MAGIC",
                                "DeepImpute", "Metformin-121"))
    


# Create UMAP plot
p <- ggplot(ds3_umap_data, aes(x = UMAP_1, y = UMAP_2, color = scale(nCount_RNA))) +
  geom_point(size = 0.3) +
  viridis::scale_color_viridis() +
  facet_grid(proportion ~ team) +
  theme_bw(base_size = 12) +
  labs(title = "UMAP Signal Correction - DS3 reads only", x = "UMAP 1", y = "UMAP 2", color = "Scaled Library Size") +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 14, face = "italic"))

ggsave("UMAP_signal_correction.png", plot = last_plot(), width = 20, height = 8)


# Upload to synapse------------------------------------------------------------------
source("utils/setup.R")
file <- synapseclient$File("UMAP_signal_correction.png", parent = "syn52812201")
stored <- syn$store(file)
