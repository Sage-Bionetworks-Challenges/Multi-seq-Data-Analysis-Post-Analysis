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
          overall_rank,
          submission_runtime,
          submission_max_memory
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
    
    if (task_n == 1) {
      # add post submissions
      post_query <- stringr::str_glue(
        "
        SELECT
          re_id AS \"id\",
          re_submitterid AS \"submitterid\",
          re_evaluationid AS \"evaluationid\",
        	prediction_fileid,
          submission_scores,
        	primary_metric_average, 
          secondary_metric_average,
          primary_rank,
          secondary_rank,
          overall_rank,
          re_submission_phase AS \"submission_phase\",
        	submission_runtime,
          submission_max_memory
        FROM {view_id}
        WHERE submission_status = 'SCORED'
          AND status = 'ACCEPTED'
          AND re_id is NOT NULL
        "
      )
      
      post_sub_df <- get_ranked_submissions(syn, post_query) %>% 
        filter(id == "9736569" | submission_phase == "public" & evaluationid == "9615324" ) %>%
        filter(team %in% c("zoradeng", "BBKCS", "moqri", "USF biostat")) %>%
        mutate(team = paste0(team, "_post")) %>%
        select(-submission_phase)
      
      sub_df <- rbind(sub_df, post_sub_df)
    }
    
    
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
