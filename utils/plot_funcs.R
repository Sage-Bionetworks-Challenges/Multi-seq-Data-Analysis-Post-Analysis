
# Boxplot of Bayes Factor on boostrapped results --------------------------
# bootstrap_boxplot <- function(.data,
#                               .x,
#                               .y,
#                               bf_col, # bayes factor column to evaluate
#                               bf_cutoffs=c(3),
#                               ref_label = "reference"
#                               ) {
#   
#   .x <- enquo(.x)
#   .y <- enquo(.y)
#   bf_col <- enquo(bf_col)
#   if (length(bf_cutoffs) > 2) stop("Too many cutoffs")
#   stopifnot(any(bf_cutoffs >= 0))
#   cutoff1 <- min(bf_cutoffs)
#   cutoff2 <- max(bf_cutoffs)
#   
#   
#   if (cutoff1 != cutoff2) {
#     r1 <- paste0("< ", cutoff1)
#     r2 <- paste0(cutoff1, " - ", cutoff2)
#     r3 <- paste0("> ", cutoff2)
#     
#     .data <- .data %>% 
#       mutate(groups = factor(
#         case_when(
#           !!bf_col == 1 ~ ref_label,
#           !!bf_col > 1 & !!bf_col <= cutoff1 ~ r1,
#           !!bf_col > cutoff1 & !!bf_col <= cutoff2 ~ r2,
#           !!bf_col > cutoff2 ~ r3,
#           TRUE ~ NA
#         ), 
#         levels = c(ref_label, r1, r2, r3)))
#     col_values <- c("#A81A50", '#F94551', "#FCB335", "#32A0B5")
#     names(col_values) <- c(ref_label, r1, r2, r3)
#   } else {
#     r1 <- paste0("< ", cutoff1)
#     r2 <- paste0("> ", cutoff1)
#     .data <- .data %>% 
#       mutate(groups = factor(
#         case_when(
#           !!bf_col == 1 ~ ref_label,
#           !!bf_col > 1 & !!bf_col <= cutoff1 ~ r1,
#           !!bf_col > cutoff1 ~ r2,
#           TRUE ~ NA
#         ), 
#         levels = c(ref_label, r1, r2)))
#     col_values <- c("#A81A50", '#F94551', "#32A0B5")
#     names(col_values) <- c(ref_label, r1, r2)
#   }
#   
#   .data %>% 
#     ggplot(aes(!!.x, !!.y, color = groups)) + 
#     labs(x = quo_name(.x), y = quo_name(.y), color = "Bayes factor") +
#     geom_boxplot(lwd = 1.2, fatten = 1) + 
#     scale_x_discrete(limits=rev) +
#     coord_flip() +
#     theme_bw(base_size = 16) + 
#     scale_color_manual(values = col_values, drop = FALSE)
# }
bootstrap_boxplot <- function(.data,
                              .x,
                              .y,
                              bf_col, # bayes factor column to evaluate
                              bf_cutoffs=c(3),
                              ref,
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
    r1 <- paste0("BF < ", cutoff1)
    r2 <- paste0(cutoff1, " < BF < ", cutoff2)
    r3 <- paste0("BF > ", cutoff2)
    
    .data <- .data %>% 
      mutate(groups = factor(
        case_when(
          !!bf_col > 0 & !!bf_col <= cutoff1 ~ r1,
          !!bf_col > cutoff1 & !!bf_col <= cutoff2 ~ r2,
          !!bf_col > cutoff2 ~ r3,
          TRUE ~ NA
        ), 
        levels = c(r1, r2, r3)))
    col_values <- c("#F94551", "#FCB335", "grey")
    names(col_values) <- c(r1, r2, r3)
  } else {
    r1 <- paste0("BF < ", cutoff1)
    r2 <- paste0("BF > ", cutoff1)
    .data <- .data %>% 
      mutate(groups = factor(
        case_when(
          !!bf_col > 0 & !!bf_col <= cutoff1 ~ r1,
          TRUE ~ NA
        ), 
        levels = c(r1, r2)))
    col_values <- c('#F94551', "grey")
    names(col_values) <- c(r1, r2)
  }
  
  ref_label_y <- .data %>% filter(!!.x == ref) %>% pull(!!.y) %>% max()
  ref_label_data <- .data %>% filter(!!.x == ref & !!.y == ref_label_y)
  .data %>% 
    ggplot(aes(!!.x, !!.y, color = groups)) + 
    labs(x = quo_name(.x), y = quo_name(.y), color = NULL) +
    geom_boxplot(lwd = 1.2, fatten = 1) + 
    scale_x_discrete(limits=rev) +
    scale_color_manual(values = col_values, drop = FALSE) +
    geom_text(data = ref_label_data, x = ref, y = ref_label_y, label = ref_label,
              hjust = -0.4, vjust = 0.5, size = 5) +
    # remove legend for geom_text
    guides(color = guide_legend(override.aes = list(label = NULL, size = NULL))) +
    coord_flip() +
    theme_bw(base_size = 16)

}

# Line plot of Ensembled models --------------------------
ensemble_ranks_line <- function(.data, 
                                .x, 
                                .y, 
                                .group, # group column used to group lines
                                linetype=TRUE
                                ) {
  
  .x <- enquo(.x)
  .y <- enquo(.y)
  .group <- enquo(.group)

  # get group that is at peak of lines in 1st group assuming it's primary group
  all_groups <- levels(as.factor(.data[[quo_name(.group)]]))
  max_x <- .data %>%
    filter(!!.group == all_groups[1]) %>%
    slice_max(order_by = !!.y) %>%
    pull(!!.x)

  # filter the data that will be plotted with points
  all_x <- levels(as.factor(.data[[quo_name(.x)]]))
  point_data <- .data %>%
    filter(!!.x %in% all_x[1:(which(all_x == max_x))])
  
  if (linetype) {
    p <- ggplot(.data, aes(!!.x, !!.y, group = !!.group, linetype = !!.group)) +
      geom_line(linewidth = 1, color = "black")
  } else {
    p <- ggplot(.data, aes(!!.x, !!.y, group = !!.group)) +
      geom_line(linewidth = 1)
  }
  p +
    geom_point(data = point_data, aes(color = !!.x), size = 6, alpha = 1) +
    scale_x_discrete(guide = guide_axis(angle = 45)) + 
    labs(x = str_to_title(quo_name(.x)), y = str_to_title(quo_name(.y)), color = str_to_title(quo_name(.x)), linewidth = str_to_title(quo_name(.group))) + 
    theme_bw(base_size = 16) 
}

