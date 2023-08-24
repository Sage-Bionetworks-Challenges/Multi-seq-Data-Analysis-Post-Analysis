# Install all R packages if not yet --------------------------
source("install.R")

# Set up synapse ----------------------------------------------------------
reticulate::use_condaenv('synapse')
synapseclient <- reticulate::import('synapseclient')
syn <- synapseclient$Synapse()
syn$login(silent = TRUE)


# Set up task variables --------------------------------------------------------
view_id <- "syn51157023"
eval_id <- list("9615023", "9615024")
gs_id <- list("syn34612394", "syn35294386")
metrics_lookup <- list(c("nrmse_score", "spearman_score"), 
                       c("summed_score", "jaccard_similarity"))


# Set up cores for parallization ------------------------------------------
ncores <- parallel::detectCores() - 1
message("<<<< ", ncores, " cores will be used for parallel computing if applicable", " <<<<")



# Set up output directory -------------------------------------------------
data_dir <- "data"
dir.create(data_dir, showWarnings = FALSE)
temp_dir <- "temp"
dir.create(temp_dir, showWarnings = FALSE)

