## Function to calculate scores on scRNAseq data ------------------
.prepare <- function(true, pred, pseudobulk = FALSE, proportion = NULL) {
  shared_cells <- intersect(colnames(pred), colnames(true))
  shared_genes <- intersect(rownames(pred), rownames(true))
  
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

calculate_nrmse <- function(.data, pseudobulk = FALSE) {
  true <- .data$true
  pred <- .data$pred
  n_na_genes <- .data$n_na_genes
  n_na_cells <- .data$n_na_cells
  total_cells <- .data$total_cells
  
  if (pseudobulk) {
    true_rs <- rowSums(true)
    pred_rs <- rowSums(pred)
    rmse <- sqrt(mean((true_rs - pred_rs)**2))
    if (n_na_genes > 0) rmse <- c(rmse, rep(1, n_na_genes))
    if (n_na_cells > 0) rmse <- rmse / (1 - n_na_cells / total_cells)
    nrmse <- rmse / (max(true_rs) - min(true_rs))
  } else {
    rmse <- sqrt(rowMeans((true - pred)**2))
    if (n_na_genes > 0) rmse <- c(rmse, rep(1, n_na_genes))
    if (n_na_cells > 0) rmse <- rmse / (1 - n_na_cells / total_cells)
    range_rr <- matrixStats::rowMaxs(true) - matrixStats::rowMins(true)
    nrmse <- rmse / range_rr
  }
  
  score <- mean(nrmse)
  
  return(score)
}

calculate_spearman <- function(.data, pseudobulk = FALSE, na.rm = TRUE) {
  true <- .data$true
  pred <- .data$pred
  n_na_genes <- .data$n_na_genes
  n_na_cells <- .data$n_na_cells
  total_cells <- .data$total_cells
  
  if (pseudobulk) {
    true_rs <- rowSums(true)
    pred_rs <- rowSums(pred)
    spearman <- cor.test(true_rs, pred_rs, method = "spearman")$estimate
    if (n_na_genes > 0) spearman <- c(spearman, rep(0, n_na_genes))
  } else {
    spearman <- sapply(1:nrow(true), function(i) cor.test(true[i, ], pred[i, ], method = "spearman")$estimate)
    if (n_na_genes > 0) spearman <- c(spearman, rep(0, n_na_genes))
  }
  
  score <- mean(spearman, na.rm = na.rm)
  if (n_na_cells > 0) score <- score * (1 - n_na_cells / total_cells)
  
  return(score)
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