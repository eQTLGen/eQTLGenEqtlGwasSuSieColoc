#!/usr/bin/env Rscript


# Load libraries
library(argparse)
library(data.table)
library(tidyverse)
library(arrow)

# Declare constants
setDTthreads(1)
parser <- ArgumentParser(description = "Filter coloc results.")

parser$add_argument("--coloc",
                    metavar = "file", type = "character",
                    help = "Coloc files."
)

parser$add_argument("--reference",
                    metavar = "file", type = "character",
                    help = "eQTLGen SNP reference file in parquet format."
)

parser$add_argument("--ld_folder",
                    metavar = "file",
                    help = "Folder containing permuted LD files."
)

# Declare function definitions

# Extract locus as data table
get_ld_matrix_wide <- function(permuted_dataset, variants) {
  start.time <- Sys.time()

  rho_mat <- permuted_dataset %>%
    filter(variant_index %in% variants) %>%
    collect() %>%
    as.data.table() %>%
    as.matrix(rownames = 1)

  start.time <- Sys.time()

  rho_mat <- rho_mat - rowMeans(rho_mat)
  # Standardize each variable
  rho_mat <- rho_mat / sqrt(rowSums(rho_mat^2))
  # Calculate correlations
  ld_matrix <- tcrossprod(rho_mat)

  end.time <- Sys.time()
  time.taken <- end.time - start.time

  print(time.taken)
  return(ld_matrix)
}


calculate_ld_for_coloc <- function(dat, variant_reference, ld_panel) {

  coloc_vars <- unique(c(dat$hit1, dat$hit2))

  coloc_variants <- variant_reference %>%
    filter(variant %in% coloc_vars) %>% collect()

  focal_chromosome <- as.integer(dat$chromosome[1])

  ld_mat <- get_ld_matrix_wide(
    ld_panel %>% filter(chr == focal_chromosome),
    coloc_variants %>% pull(variant_index))

  variant_map <- deframe(coloc_variants %>% select(variant, variant_index))

  dat$ld_R <- ld_mat[variant_map[dat$hit1], variant_map[dat$hit2]]

  return(dat)
}


# Main

#' Execute main
#' 
#' @param argv A vector of arguments normally supplied via command-line.
main <- function(argv = NULL) {
  if (is.null(argv)) {
    argv <- commandArgs(trailingOnly = T)
  }

  args <- parser$parse_args(argv)

  # Process input
  coloc_results <- fread(args$coloc) %>%
    separate_wider_delim(eQTL_locus, delim=":", names=c("chromosome", "basepairs"), cols_remove=F)

  variant_reference <- open_dataset(args$reference)

  # Perform method
  ld_folder <- args$ld_folder
  ld_panel <- open_dataset(ld_folder)

  coloc_with_ld <- coloc_results %>%
    as.data.table() %>%
    group_by(chromosome) %>%
    group_modify(~ calculate_ld_for_coloc(.x, variant_reference, ld_panel), .keep = T)

  fwrite(coloc_with_ld, "coloc_results_with_ld_col.txt", sep="\t", row.names = F, col.names=T, quote=F)
  # Process output
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}