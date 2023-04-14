suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
})

message("\n")
message("------------------------------------------")
message("The script is only working for scRNA-seq and scATAC-seq Data Analysis DREAM Challenge")
message("------------------------------------------")

source("utils.R")

# set up synapse
reticulate::use_condaenv('synapse')
synapseclient <- reticulate::import('synapseclient')
challengeutils <- reticulate::import("challengeutils")
syn <- synapseclient$Synapse()
syn$login(silent = TRUE)

# set args
parser <- argparse::ArgumentParser()
parser$add_argument("--task", help = "task1 or task2")
parser$add_argument("--phase", help = "public or private")
parser$add_argument("--dry-run", action="store_true", help = "on/off dryrun")
args <- parser$parse_args()


view_id <- "syn51157023"
tasks <- list(task1 = list(old_eval = "9615023", new_eval = "9615324"),
              task2 = list(old_eval = "9615024", new_eval = "9615302"))
admin_names <- c("scRNA-seq and scATAC-seq Data Analysis DREAM Challenge Organizers", "rchai")

task <- args$task
stopifnot(task %in% c("task1", "task2"))
phase <- args$phase
stopifnot(phase %in% c("public", "private"))
dry_run <- args$dry_run

# get ids of all scored submissions
message("Downloading submissions ...")
eval_id <- tasks[[task]]$old_eval
task_sub_df <- get_ranked_submissions(syn, view_id, eval_id, phase)

# validate if submitter is overlapped with existing teams
all_users <- unique(task_sub_df$submitterid)
all_teams <- all_users[sapply(all_users, is_team, syn = syn)]
is_dup_user <- validate_users(syn,
                              users = all_users,
                              teams = all_teams,
                              .drop = TRUE)

if (length(is_dup_user) > 0) {
  message(
    "Submissions from submitters will not be re-submitted:\n",
    paste(task_sub_df$team[match(names(is_dup_user), task_sub_df$submitterid)], collapse = '", "'),
    "\nbecause of duplicated with exiting teams"
  )
}

# remove duplicated user and get the best submission of a team/user
task_sub_df <- task_sub_df %>%
  filter(!submitterid %in% names(is_dup_user),
         !team %in% admin_names) %>%
  mutate(overall_rank = as.numeric(overall_rank)) %>%
  group_by(team) %>%
  slice_min(overall_rank, n = 1) %>%
  arrange(overall_rank)

if (nrow(task_sub_df) > 0) {

  message("\nRe-submitting submissions from evluation ", dQuote(eval_id), " and ", phase, " phase ...")

  sub_ids <- unique(task_sub_df$id)
  for (n in seq_along(sub_ids)) {

    message(str_glue("Submitting {dQuote(sub_ids[n])} ({n} of {nrow(task_sub_df)}) ..."))

    if (!dry_run) {

      # this approach not works on the project lack of access
      # submitted <- resubmit(syn, sub_ids[n], tasks[[task]]$new_eval)

      # copy images to staging then submit
      staging_project_id <- "syn26720921"
      sub <- syn$getSubmission(sub_ids[n])
      sub_image <- str_glue("{sub$dockerRepositoryName}@{sub$dockerDigest}")
      image_tag <- str_glue("{sub$id}_{phase}")

      image <- copy_model(image = sub_image,
                          project_id = staging_project_id,
                          name = "post-submissions",
                          tag = image_tag)

      repo_entity <- syn$restGET(str_glue('/entity/dockerRepo/id?repositoryName={image$repo_name}'))

      submitted <- syn$submit(
        evaluation = tasks[[task]]$new_eval,
        entity = repo_entity$id,
        name = sub$id,
        dockerTag = image_tag
      )

      # wait 15 min to check if submission is scored - not working
      # cp <- 0
      # while (cp <= 0) {
      #   Sys.sleep(15 * 60)
      #   status <- syn$getSubmissionStatus(submitted$id)
      #   is_score <- tryCatch(
      #     {
      #       res <- status$submissionAnnotations$submission_status == "SCORED" ||
      #         status$status == "INVALID"
      #       return(res)
      #     },
      #     error = function(e) FALSE
      #   )
      #   if (is_score) cp <- cp + 1
      # }
      
      # update annots with origin submitterid and phase
      Sys.sleep(60) # wait 1 min to update annotations
      message("Annotating original submitterid and submission phase ... ")
      # id and submitterid cannot be overwrote, so add new columns
      team_id <- task_sub_df$submitterid[match(sub_ids[n], task_sub_df$id)]
      annots <- list(
        re_submission_phase = as.character(phase), # original submission_phase
        re_submitterid = as.character(team_id) # original submitterid
      )
      annot_res <- challengeutils$annotations$annotate_submission(syn, submitted$id, annots)
    }
  }
}
