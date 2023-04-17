# Check if a synapse team --------------------------
is_team <- function(syn, uid) {
  if (missing(syn)) stop('argument "syn" is missing')
  if (missing(uid)) stop('argument "uid" is missing')
  out <- tryCatch(
    {
      syn$getTeam(uid)
      return(TRUE)
    },
    error = function(e) FALSE
  )
  return(out)
}


# Check if a synapse user --------------------------
is_user <- function(syn, uid) {
  if (missing(syn)) stop('argument "syn" is missing')
  if (missing(uid)) stop('argument "uid" is missing')
  out <- tryCatch(
    {
      syn$getUserProfile(uid)
      return(TRUE)
    },
    error = function(e) FALSE
  )
  return(out)
}


# Check if a user is part of an existing team --------------------------
is_member <- function(syn, user_uid, team_uid) {
  if (missing(user_uid)) stop('argument "user_uid" is missing')
  if (missing(team_uid)) stop('argument "team_uid" is missing')
  tryCatch(
    {
      member_uids <- reticulate::iterate(syn$getTeamMembers(team_uid)) %>%
        sapply(., function(member_object) {
          member <- jsonlite::fromJSON(as.character(member_object))
          return(member$member$ownerId)
        })
      if (length(member_uids) == 0) {
        return(FALSE)
      } else {
        return(user_uid %in% member_uids)
      }
    },
    error = function(e) FALSE
  )
}


# Check if users are part of existing teams --------------------------
validate_users <- function(syn, users, teams, .drop = FALSE) {
  res <- sapply(users, function(user) {
    if (is_user(syn, user)) { # only validate user
      
      team_res <- sapply(teams, function(team) {
        is_member(syn = syn, user_uid = user, team_uid = team)
      })
      return(any(team_res))
    } else {
      return(FALSE)
    }
  })
  
  if (.drop) {
    res <- res[res]
  }
  return(res)
}


# Retrieve the user/team display name --------------------------
get_name <- function(syn, id) {
  name <- tryCatch(
    {
      syn$getUserProfile(id)$userName
    },
    error = function(err) {
      syn$getTeam(id)$name
    }
  )
  return(name)
}


# Resubmit docker models using synapseclient --------------------------
# Not working for the projects without access
# resubmit <- function(syn, sub_id, new_eval_id) {
#   # resubmit the model without minimal requirements for the request
#   # no teamId, contributors or eligibilityStateHash
#   # note, it only works for the projects with write access
#   sub <- syn$getSubmission(sub_id)
#
#   submission <- list(
#     'evaluationId' = new_eval_id,
#     'name' = sub_id,
#     'entityId' = sub["entityId"],
#     'versionNumber' = sub$get('versionNumber', 1),
#     'dockerDigest' = sub["dockerDigest"],
#     'dockerRepositoryName' = sub["dockerRepositoryName"],
#     'teamId' = NULL,
#     'contributors' = NULL,
#     'submitterAlias' = NULL
#   )
#
#   docker_repo_entity <- syn$restGET(stringr::str_glue('/entity/dockerRepo/id?repositoryName={sub["dockerRepositoryName"]}'))
#   entity <- syn$get(docker_repo_entity["id"], downloadFile=FALSE)
#   uri <- stringr::str_glue("/evaluation/submission?etag={entity['etag']}")
#
#   # ignore eligibility
#   # eligibility <- syn$restGET(stringr::str_glue('/evaluation/{sub["evaluationId"]}/team/{sub["teamId"]}/submissionEligibility'))
#   # uri <- stringr::str_glue("{uri}&submissionEligibilityHash={eligibility['eligibilityStateHash']}")
#   submitted <- syn$restPOST(uri, jsonlite::toJSON(submission, auto_unbox = TRUE, null = "null"))
#   return(submitted)
# }


# Copy and collect docker models to other project --------------------------
copy_model <- function(image, project_id, name, tag = "latest") {
  
  # get new project repo
  docker_repo <- stringr::str_glue("docker.synapse.org/{project_id}")
  
  # TODO: add validation on image string
  
  # get docker image names
  repo_name <- file.path(docker_repo, name)
  new_image <- stringr::str_glue("{repo_name}:{tag}")
  
  system(stringr::str_glue("docker pull {image}"))
  system(stringr::str_glue("docker tag {image} {new_image}"))
  system(stringr::str_glue("docker push {new_image}"))
  system(stringr::str_glue("docker image rm {image} {new_image}"))
  
  return(list(repo_name = repo_name, tag = tag))
}


get_ranked_submissions <- function(syn, query) {

  # download the submissions ordered by overall rank
  sub_df <- syn$tableQuery(query)$asDataFrame() %>%
    mutate(across(everything(), as.character),
           team = as.character(sapply(submitterid, get_name, syn = syn)))
  return(sub_df)
}


# Retrieving scores from submission view table --------------------------
get_scores <- function(syn, sub_df) {
  # validate if any valid submission to prevent from failing
  stopifnot(nrow(sub_df) > 0)
  stopifnot("submission_scores" %in% colnames(sub_df))
  
  # read all valid scores results
  all_scores <- lapply(1:nrow(sub_df), function(sub_n) {
    score_id <- sub_df$submission_scores[sub_n]
    
    # read all test case scores for each submission
    score_df <- syn$get(score_id)$path %>%
      data.table::fread(data.table = FALSE, verbose = FALSE) %>%
      mutate(
        id = sub_df$id[sub_n],
        submitterid = sub_df$submitterid[sub_n],
        team = sub_df$team[sub_n]
      )
    return(score_df)
  }) %>% bind_rows()
}


# Rank metrics across all submissions --------------------------
rank_submissions <- function(scores, primary_metric, secondary_metric, group=c("id", "team")) {
  stopifnot(nrow(scores) > 0)
  stopifnot(c(primary_metric, secondary_metric, "dataset") %in% colnames(scores))
  stopifnot(group %in% colnames(scores))
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
    group_by_at(group) %>%
    # get average scores of all testcases ranks in one submission
    summarise(
      primary_rank = mean(testcase_primary_rank),
      secondary_rank = mean(testcase_secondary_rank),
      .groups = "drop"
    ) %>%
    # rank overall rank on primary, tie breaks by secondary
    arrange(primary_rank, secondary_rank) %>%
    mutate(overall_rank = row_number())
  
  return(rank_df)
}
