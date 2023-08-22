# Bootstrap on rows of matrix and return a long format result --------------------------
simple_bootstrap <- function(.data,
                             seq_size,
                             .by = NULL,
                             n_iter = 1000,
                             seed = 98109,
                             ncores = 1) {
  
  # generate the re-sampling indices
  set.seed(seed)
  bootstrap_indices <- matrix(sequence(seq_size), seq_size, n_iter) %>%
    apply(2, sample, replace = TRUE)
  
  bs_results <- parallel::mclapply(1:ncol(bootstrap_indices), function(n) {
    iter_indice <- bootstrap_indices[, n]
    rs_data <- .data %>%
      slice(iter_indice, .by = !!sym(.by)) %>%
      mutate(bs_n = n)
    
    return(rs_data)
  }, mc.cores = ncores) %>% bind_rows()
  
  return(bs_results)
}


# Bootstrap data matrix and return a wide format result --------------------------
# A custom function can be applied to each bootstrapped result
bootstrap <- function(.data,
                      func = NULL,
                      ...,
                      direction = 1,
                      n_iter = 1000,
                      seed = 98109,
                      ncpu = 1) {

  stopifnot(is.numeric(direction) & direction %in% c(1, 2))

  seq_size <- ifelse(direction == 1, nrow(.data), ncol(.data))
  
  # generate the re-sampling indices
  set.seed(seed)
  bootstrap_indices <- matrix(sequence(seq_size), seq_size, n_iter) %>%
    apply(2, sample, replace = TRUE)

  pbmcapply::pbmclapply(1:ncol(bootstrap_indices), function(i) {
    
    # resample
    iter_indice <- bootstrap_indices[, i]

    if (direction == 1) {
      bs_data <- .data[iter_indice, ]
    } else {
      bs_data <- .data[, iter_indice]
    }
    
    # apply the custom function on bootstrapped data,
    # otherwise, just return bootstrapped data
    if (!is.null(func)) {
      res <- func(bs_data, ...)
    } else {
      res <- bs_data
    }

    if (!is.atomic(res)) stop("Error: The output of ",
                              deparse(substitute(func)),
                              " must be one-dimensional vector.")
    return(res)
  }, mc.cores = ncpu, ignore.interactive = TRUE) %>%
    setDT(unlist(., recursive = FALSE)) %>%
    data.table::setnames(paste0("bs_", sequence(n_iter)))
}


# Compute Bayes Factor between two vectors --------------------------
bayes_factor <- function(model, ref) {
  
  # if same values, return 0
  if (identical(model, ref)) return(0)
  
  diff <- model - ref
  
  pos_n <- sum(diff >= 0)
  neg_n <- sum(diff < 0)
  
  # prevent from all zeros
  K <- (pos_n + 0.01) / (neg_n + 0.01)
  
  # reciprocate fraction of K
  if (K < 1) K <- 1 / K
  
  return(K)
}


