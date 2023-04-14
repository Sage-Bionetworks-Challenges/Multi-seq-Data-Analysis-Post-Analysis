# instead of using renv, use pacman
# since I got issue when restore Seurat using renv
# and not many packages need to be installed,

if (!require("pacman", quietly = TRUE)) install.packages("pacman")

pkgs <- c("dplyr", "stringr", 
          "data.table", 
          "ggplot2", "patchwork", 
          "Seurat", "bedr")

pacman::p_load(pkgs, character.only = TRUE)

