get_name <- function(id) {
  name <- tryCatch({
    syn$getUserProfile(id)$userName
  }, error = function(err) {
    syn$getTeam(id)$name
  })
  name
}


get_submissions <- function(view_id, gs_id) {
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
      prediction_fileid
    FROM {sub_id} 
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
           team = sapply(submitterid, get_name),
    )
  sub_df
}


bootstrap <- function(seq_size,
                      n_iterations=1000,
                      seed=98109) {
  set.seed(seed)
  bs_indices <- lapply(sequence(n_iterations), function(n) {
    sample(sequence(seq_size), seq_size, replace = TRUE)
  })
  bs_indices
}


rank_submission <- function(scores, # matrix of scores for all testcases all submissions
                            dataset_col, # column containing the test case file names 
                            primary_col, # primary metric column name
                            bs_col, # column containing the bootstrapping indice
                            secondary_col # secondary metric column name
                            ) {
  # should contain the testcase
  stopifnot(c("id", "submitterid", dataset_col, bs_col, primary_col) %in% colnames(scores))
  message("Ranking scores ...")
  # rank the scores
  rank_df <-
    scores %>%
    group_by(all_of(dataset_col, bs_col)) %>%
    # rank each testcase score of one submission compared to all submissions
    # the smaller values, the smaller ranks, aka higher ranks
    mutate(
      testcase_primary_rank = rank(-!!primary_col),
      testcase_secondary_rank = rank(-!!secondary_col)
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
  rank_df
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