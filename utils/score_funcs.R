## Function to calculate scores on scRNAseq data ------------------
.prepare <- function(true, pred, pseudobulk = FALSE, proportion = NULL) {
  shared_cells <- intersect(colnames(pred), colnames(true))
  shared_genes <- intersect(rownames(pred), rownames(true))

  if (is.data.frame(true)) true <- as.matrix(true)
  # number of missing cells and genes
  if (pseudobulk) {
    stopifnot(!is.null(proportion))
    proportion <- as.numeric(gsub("p0", "0.", proportion)) # i.e, 'p025' -> 0.25
    # using pseudobulk means the cells of training data are subsetted
    # calculate expected number of cells
    total_cells <- floor(ncol(true) * proportion)
    n_na_cells <- total_cells - length(shared_cells)
  } else {
    n_na_cells <- sum(!colnames(true) %in% colnames(pred))
    total_cells <- ncol(true)
  }
  n_na_genes <- sum(!rownames(true) %in% rownames(pred))
  total_genes <- nrow(true)
  
  # match the orders of genes and cells
  if (pseudobulk) {
    out_true <- as.matrix(true[shared_genes, ])
    out_pred <- as.matrix(pred[shared_genes, ])
  } else {
    out_true <- as.matrix(true[shared_genes, shared_cells])
    out_pred <- as.matrix(pred[shared_genes, shared_cells])
  }
  return(list(
    true = out_true,
    pred = out_pred,
    n_na_cells = n_na_cells,
    n_na_genes = n_na_genes,
    total_cells = total_cells,
    total_genes = total_genes
  ))
}


calculate_nrmse <- function(.data, pseudobulk = FALSE, penality = TRUE, aggregate_func = NULL, ...) {
  true <- .data$true
  pred <- .data$pred
  n_na_genes <- .data$n_na_genes
  n_na_cells <- .data$n_na_cells
  total_cells <- .data$total_cells
  
  if (pseudobulk) {
    true_rs <- rowSums(true)
    pred_rs <- rowSums(pred)
    rmse <- sqrt(mean((true_rs - pred_rs)**2))
    if (n_na_genes > 0 & penality) rmse <- c(rmse, rep(1, n_na_genes))
    if (n_na_cells > 0 & penality) rmse <- rmse / (1 - n_na_cells / total_cells)
    nrmse <- rmse / (max(true_rs) - min(true_rs))
  } else {
    rmse <- sqrt(rowMeans((true - pred)**2))
    if (n_na_genes > 0 & penality) rmse <- c(rmse, rep(1, n_na_genes))
    if (n_na_cells > 0 & penality) rmse <- rmse / (1 - n_na_cells / total_cells)
    range_rr <- matrixStats::rowMaxs(true) - matrixStats::rowMins(true)
    nrmse <- rmse / range_rr
    nrmse[which(is.infinite(nrmse))] <- NaN
  }
  
  if (!is.null(aggregate_func)) {
    score <- aggregate_func(nrmse, ...)
  } else {
    score <- nrmse
  }

  return(score)
}


calculate_spearman <- function(.data, pseudobulk = FALSE, penality = TRUE, aggregate_func = NULL, ..., ncores = 1) {
  true <- .data$true
  pred <- .data$pred
  n_na_genes <- .data$n_na_genes
  n_na_cells <- .data$n_na_cells
  total_cells <- .data$total_cells
  
  if (pseudobulk) {
    true_rs <- rowSums(true)
    pred_rs <- rowSums(pred)
    spearman <- cor.test(true_rs, pred_rs, method = "spearman")$estimate
    if (n_na_genes > 0 & penality) spearman <- c(spearman, rep(0, n_na_genes))
  } else {
    spearman <- simplify2array(
      parallel::mclapply(1:nrow(true), mc.cores = ncores, function(i) cor.test(true[i, ], pred[i, ], method = "spearman")$estimate)

    )
    if (n_na_genes > 0 & penality) spearman <- c(spearman, rep(0, n_na_genes))
  }
  
  if (!is.null(aggregate_func)) {
    score <- aggregate_func(spearman, ...)
  } else {
    score <- spearman
  }

  if (n_na_cells > 0 & penality) score <- score * (1 - n_na_cells / total_cells)
  return(as.numeric(score))
}


## Function to calculate scores on scATACseq data ------------------
category_recall <- function(a, b, verbose = FALSE) {
  a.int2 <- tryCatch(
    {
      bedr(input = list(a = a, b = b), method = "intersect -u", params = "-sorted", verbose = verbose)
    },
    error = function(cond) {
      return(data.frame())
    } # return empty df on error
  )
  # fraction of rows in ground truth that overlap at all with peaks in submission
  return(nrow(a.int2) / nrow(a))
}


## Wrap up funcion to scoring bootstrapped scRNAseq data imputation  ------------------
score_bs_imputation <- function(bs_imputed, # bootstrapped imputed
                                gs, # groundtruth
                                proportion,
                                n_na_cells,
                                n_na_genes, 
                                total_cells,
                                total_genes,
                                pseudobulk = FALSE) {
  
  # match the order of genes and cells of grouth truth with bootraped data
  stopifnot(any(rownames(bs_imputed) %in% rownames(gs), colnames(bs_imputed) %in% colnames(gs)))
  bs_gs <- gs[rownames(bs_imputed), colnames(bs_imputed)]

  # not re-normalize the boostrapped data since it will change the origin variability
  
  # prepare the data for scoring
  # set the numbers from origin data (not bootstrapped data) for penality
  # aka, perserve the penality of origin imputed data if any missing genes/cells
  bs_eval_data <- list(
    true = bs_gs,
    pred = bs_imputed,
    n_na_cells = n_na_cells,
    n_na_genes = n_na_genes,
    total_cells = total_cells,
    total_genes = total_genes
  )
  # scoring
  nrmse_score <- calculate_nrmse(bs_eval_data, pseudobulk = pseudobulk, aggregate_func = mean, na.rm = TRUE)
  spearman_score <- calculate_spearman(bs_eval_data, pseudobulk = pseudobulk, aggregate_func = mean, na.rm = TRUE)

  return(as.numeric(c(nrmse_score, spearman_score)))
}


## Function to wrap scoring functions for scATACseq data ------------------
score_bs_peaks <- function(bs_peak, # bootstrapped call peaks
                           gs, # grountruth data for evaluation 
                           use_cellspecific = FALSE) {
  
  # retrieve peaks from groundtruth data will be used to score
  gs_ranked_filtered <- gs[["gs_ranked_filtered"]]
  gs_sort <- gs[["gs_sort"]]
  gs_sort_ubi <- gs[["gs_sort_ubi"]]
  gs_sort_tss <- gs[["gs_sort_tss"]]
  gs_sort_csp <- gs[["gs_sort_csp"]]

  # sort bootstrapped peaks
  peak_sort <- bedr.merge.region(
    bedr.sort.region(
      as.data.frame(bs_peak),
      check.zero.based = FALSE, 
      check.chr = FALSE, 
      check.valid = FALSE,
      verbose = FALSE
    ), 
    verbose = FALSE
  )

  # caculate scores
  j_result <- jaccard(gs_sort, peak_sort, check.merge = TRUE, check.chr = FALSE, check.sort = FALSE, check.valid = FALSE, verbose = FALSE)
  jaccard_similarity <- as.numeric(j_result$jaccard)
  recall_ubiquitous <- category_recall(gs_sort_ubi, pred_sort)
  recall_tss <- category_recall(gs_sort_tss, pred_sort)
  recall_cellspecific <- category_recall(gs_sort_csp, pred_sort)
  
  summed_score <- jaccard_similarity + recall_ubiquitous + recall_tss
  
  # report cellspecific for ds2, but not added to summed score
  if (use_cellspecific) summed_score <- summed_score + recall_cellspecific

  return(c(
    jaccard_similarity, 
    recall_ubiquitous, 
    recall_tss,
    recall_cellspecific,
    summed_score
  ))
}

      

      