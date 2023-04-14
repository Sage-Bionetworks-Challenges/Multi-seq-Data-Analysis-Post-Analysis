# return a list of bootstrapped indices
bootstrap_indices <- function(seq_size,
                              n_iterations = 1000,
                              seed = 98109) {
  set.seed(seed)
  bs_indices <- lapply(sequence(n_iterations), function(n) {
    sample(sequence(seq_size), seq_size, replace = TRUE)
  })
  bs_indices
}


# bootstrap the data matrix in long format, 
# bootstrap <- function(.data,
#                       seq_size,
#                       .by = NULL,
#                       n_iterations = 1000,
#                       seed = 98109,
#                       ncores = 1) {
#   set.seed(seed)
  
#   rs_indices <- bootstrap_indices(seq_size, n_iterations, seed)
  
#   bs_results <- parallel::mclapply(seq_along(rs_indices), function(n) {
    
#     rs_data <- .data %>%
#       slice(rs_indices[[n]], .by = !!sym(.by)) %>%
#       mutate(bs_n = n)
    
#     return(rs_data)
#   }, mc.cores = ncores) %>% bind_rows()
  
#   return(bs_results)
# }


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


bootstrap <- function(.data,
                      func = NULL,
                      ...,
                      direction = 1,
                      n_iter = 10,
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
    
    # apply the custom function on boostrapped data,
    # otherwise, just return boostrapped data
    if (!is.null(func)) {
      res <- func(bs_data, ...)
    } else {
      res <- bs_data
    }

    if (!is.atomic(res)) stop("Error: The output of ", 
                              deparse(substitute(data)),
                              " must be one-dimensional vector.")
    return(res)
  }, mc.cores = ncpu, ignore.interactive = TRUE) %>%
    setDT(unlist(., recursive = FALSE)) %>%
    data.table::setnames(paste0("bs_", sequence(n_iter)))
}
