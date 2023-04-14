source("utils/setup.R")
source("utils/synapse_funcs.R")


# Download submissions ----------------------------------------------------
invisible(
  lapply(1:2, function(task_n) {

    # query the final submission, ordered by final rank
    query <- stringr::str_glue(
      "
        SELECT
          id,
          submitterid,
          evaluationid,
          prediction_fileid,
          submission_scores,
          primary_metric_average,
          secondary_metric_average,
          primary_rank,
          secondary_rank,
          overall_rank
        FROM {view_id}
        WHERE
          status = 'ACCEPTED'
          AND submission_status = 'SCORED'
          AND submission_phase = 'private'
          AND evaluationid = '{eval_id[[task_n]]}'
          AND overall_rank IS NOT NULL
        ORDER BY overall_rank
        "
    )
    
    sub_df <- get_ranked_submissions(syn, query)
    
    saveRDS(sub_df, 
            file.path(data_dir, str_glue("final_submissions_task{task_n}.rds")))
  })
)

# Retrieve all scores -----------------------------------------------------
# download all test case's score each submission
invisible(
  lapply(1:2, function(task_n) {
    
    sub_df <- readRDS(
      file.path(data_dir, str_glue("final_submissions_task{task_n}.rds"))
    )
    
    scores_df <- get_scores(syn, sub_df)
    
    saveRDS(scores_df, 
            file.path(data_dir, str_glue("final_scores_task{task_n}.rds"))
    )
  })
)
