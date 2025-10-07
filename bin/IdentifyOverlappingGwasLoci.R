#!/usr/bin/env Rscript


# Load libraries
library(GenomicRanges)
library(dplyr)

# Declare constants

# Declare function definitions

# Main

#' Execute main
#' 
#' @param argv A vector of arguments normally supplied via command-line.
main <- function(argv = NULL) {
  if (is.null(argv)) {
    argv <- commandArgs(trailingOnly = T)
  }

  gwas_finemap <- tibble(locus_file = list.files("GWAS_finemap"))

  locus_mat <- str_extract(gwas_finemap$locus_file, "(.+)__(\\d+):(-?\\d+)_(\\d+)___gwas.txt.gz", group=c(1,2,3,4))
  colnames(locus_mat) <- c("gwas_id", "chromosome", "start", "stop")

  locus_tibble <- as_tibble(locus_mat) %>%
    mutate(locus_file = gwas_finemap$locus_file,
           across(c(start, stop), as.integer))

  # Example dataframe
  df <- data.frame(
    chromosome = c("1", "1", "1", "2", "2"),
    start = c(100, 150, 300, 500, 520),
    stop = c(200, 250, 400, 550, 600),
    gwas_trait = c("TraitA", "TraitA", "TraitA", "TraitB", "TraitB")
  )

  # Function to get pairwise overlaps within each trait
  get_overlaps <- function(subdf) {
    gr <- GRanges(seqnames = subdf$chromosome,
                  ranges = IRanges(start = subdf$start, end = subdf$stop))

    hits <- findOverlaps(gr)  # returns all overlaps, including self

    overlaps <- data.frame(
      locus1 = queryHits(hits),
      locus2 = subjectHits(hits)
    ) %>%
      filter(locus1 < locus2) %>%   # remove self-overlaps & duplicates
      mutate(
        chr = subdf$chromosome[locus1],
        start1 = subdf$start[locus1],
        stop1  = subdf$stop[locus1],
        start2 = subdf$start[locus2],
        stop2  = subdf$stop[locus2]
      )

    return(overlaps)
  }

  # Apply per trait
  pairwise_overlaps <- locus_tibble %>%
    group_by(gwas_id) %>%
    group_modify(~ get_overlaps(.x))

  print(pairwise_overlaps)
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}