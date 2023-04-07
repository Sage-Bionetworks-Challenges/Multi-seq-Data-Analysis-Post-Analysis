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
bootstrap <- function(.data,
                      seq_size,
                      .by = NULL,
                      n_iterations = 1000,
                      seed = 98109,
                      ncores = 1) {
  set.seed(seed)
  
  rs_indices <- bootstrap_indices(seq_size, n_iterations, seed)
  
  bs_results <- parallel::mclapply(seq_along(rs_indices), function(n) {
    
    rs_data <- .data %>%
      slice(rs_indices[[n]], .by = !!sym(.by)) %>%
      mutate(bs_n = n)
    
    return(rs_data)
  }, mc.cores = ncores) %>% bind_rows()
  
  return(bs_results)
}


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