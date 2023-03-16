get_name <- function(syn, id) {
  name <- tryCatch({
    syn$getUserProfile(id)$userName
  }, error = function(err) {
    syn$getTeam(id)$name
  })
  return(name)
}


get_submissions <- function(view_id) {
  # set up synapse
  reticulate::use_condaenv('synapse')
  synapseclient <- reticulate::import('synapseclient')
  syn <- synapseclient$Synapse()
  syn$login(silent = TRUE)
  
  # query the submission view
  query <- str_glue(
    "
    SELECT 
      id,
      submitterid,
      evaluationid,
      prediction_fileid,
      submission_scores
    FROM {view_id} 
    WHERE 
      status = 'ACCEPTED'
      AND submission_status = 'SCORED' 
      AND submission_phase = 'private'
      AND overall_rank IS NOT NULL
    ORDER BY overall_rank
    "
  )
  
  # download the submissions
  sub_df <- syn$tableQuery(query)$asDataFrame() %>%
    mutate(across(everything(), as.character),
           task = if_else(evaluationid == "9615023", "task1", "task2"),
           team = sapply(submitterid, get_name, syn = syn),
    )
  return(sub_df)
}


get_scores <- function(sub_df) {
  # validate if any valid submission to prevent from failing
  stopifnot(nrow(sub_df) > 0)

  # read all valid scores results
  all_scores <- lapply(1:nrow(sub_df), function(sub_n) {

    score_id <- sub_df$submission_scores[sub_n]
    
    # read all test case scores for each submission
    score_df <- syn$get(score_id)$path %>%
      data.table::fread(verbose = FALSE) %>%
      mutate(id = sub_df$id[sub_n],
             submitterid = sub_df$submitterid[sub_n],
             team = sub_df$team[sub_n])

    return(score_df)
  }) %>% bind_rows()
}


rank_submissions <- function(scores, primary_metric, secondary_metric) {
  stopifnot(nrow(scores) > 0)
  stopifnot(c(primary_metric, secondary_metric, "dataset", "id", "submitterid") %in% colnames(scores))
  # rank the scores
  rank_df <-
    scores %>%
    group_by(dataset) %>%
    # rank each testcase score of one submission compared to all submissions
    # the smaller values, the smaller ranks, aka higher ranks
    mutate(
      testcase_primary_rank = rank(-(!!sym(primary_metric))),
      testcase_secondary_rank = rank(-(!!sym(secondary_metric)))
    ) %>%
    group_by(id, submitterid) %>%
    # get average scores of all testcases ranks in one submission
    summarise(
      avg_primary_rank = mean(testcase_primary_rank),
      avg_secondary_rank = mean(testcase_secondary_rank),
      .groups = 'drop'
    ) %>%
    # rank overall rank on primary, tie breaks by secondary
    arrange(avg_primary_rank, avg_secondary_rank) %>%
    mutate(overall_rank = row_number())
}


bootstrap <- function(.data,
                      seq_size,
                      func,
                      ...,
                      .by,
                      n_iterations=1000,
                      seed=98109,
                      ncores=1) {
  
  set.seed(seed)

  rs_indices <- lapply(sequence(n_iterations), function(n_bs) {
    sample(sequence(seq_size), seq_size, replace = TRUE)
  })

  bs_results <- parallel::mclapply(seq_along(rs_indices), function(n) {

    rs_data <- scores %>%
      slice(rs_indices[[n]], .by = !!sym(.by))

    bs_result <-  func(rs_data, ...) %>%
      mutate(n_bs = n) 

    return(bs_result)
  }, mc.cores = ncores) %>% bind_rows() 
  
  return(bs_results)
}


compute_bayes_factor <- function(bootstrapMetricMatrix,
                                 refPredIndex,
                                 invertBayes){
  
  M <- as.data.frame(bootstrapMetricMatrix - bootstrapMetricMatrix[,refPredIndex])
  K <- apply(M ,2, function(x) {
    k <- sum(x >= 0)/sum(x < 0)
    
    # Logic handles whether reference column is the best set of predictions.
    if(sum(x >= 0) > sum(x < 0)){
      return(k)
    }else{
      return(1/k)
    }
  })
  K[refPredIndex] <- 0
  if(invertBayes == T){K <- 1/K}
  return(K)
}