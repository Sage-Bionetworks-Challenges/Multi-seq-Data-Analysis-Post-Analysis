
bootstrap_boxplot <- function(.data,
                              .x,
                              .y,
                              bf_col, # bayes factor column to evaluate
                              bf_cutoffs=c(5),
                              ref_label = "reference"
                              ) {
  
  .x <- enquo(.x)
  .y <- enquo(.y)
  bf_col <- enquo(bf_col)
  if (length(bf_cutoffs) > 2) stop("Too many cutoffs")
  stopifnot(any(bf_cutoffs >= 0))
  cutoff1 <- min(bf_cutoffs)
  cutoff2 <- max(bf_cutoffs)
  
  
  if (cutoff1 != cutoff2) {
    r1 <- paste0("< ", cutoff1)
    r2 <- paste0(cutoff1, " - ", cutoff2)
    r3 <- paste0("> ", cutoff2)
    
    .data <- .data %>% 
      mutate(groups = factor(
        case_when(
          !!bf_col > 0 & !!bf_col <= cutoff1 ~ r1,
          !!bf_col > cutoff1 & !!bf_col <= cutoff2 ~ r2,
          !!bf_col > cutoff2 ~ r3,
          TRUE ~ ref_label
        ), 
        levels = c(ref_label, r1, r2, r3)))
    col_values <- c("#A81A50", '#F94551', "#FCB335", "#32A0B5")
    names(col_values) <- c(ref_label, r1, r2, r3)
  } else {
    r1 <- paste0("< ", cutoff1)
    r2 <- paste0("> ", cutoff1)
    .data <- .data %>% 
      mutate(groups = factor(
        case_when(
          !!bf_col > 0 & !!bf_col <= cutoff1 ~ r1,
          !!bf_col > cutoff1 ~ r2,
          TRUE ~ ref_label
        ), 
        levels = c(ref_label, r1, r2)))
    col_values <- c("#A81A50", '#F94551', "#32A0B5")
    names(col_values) <- c(ref_label, r1, r2)
  }
  
  .data %>% 
    ggplot(aes(!!.x, !!.y, color = groups)) + 
    labs(x = quo_name(.x), y = quo_name(.y), color = "Bayes factor") +
    geom_boxplot(lwd = 1.2, fatten = 1) + 
    scale_x_discrete(limits=rev) +
    coord_flip() +
    theme_classic(base_size = 16) + 
    scale_color_manual(values = col_values, drop = FALSE)
}



ensemble_ranks_line <- function(.data, 
                                .x, 
                                .y, 
                                .group # group column used to facet lines
                                ) {
  
  .x <- enquo(.x)
  .y <- enquo(.y)
  .group <- enquo(.group)
  
  ggplot(.data, aes(!!.x, !!.y, color = !!.group, group = !!.group)) +
    geom_line(linewidth = 1) +
    scale_x_discrete(guide = guide_axis(angle = 45)) + 
    labs(x = quo_name(.x), y = quo_name(.y), color = quo_name(.group)) +
    theme_bw(base_size = 16)
}

