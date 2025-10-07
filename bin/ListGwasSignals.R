#!/usr/bin/env Rscript


# Load libraries
library(argparse)
library(data.table)

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
  # Process input

  # Create a parser
  parser <- ArgumentParser(description = 'Read GWAS signals from a list of loci')

  # Add arguments
  parser$add_argument("--loci", required = TRUE, help = "File with finemapping loci", nargs="+")
  parser$add_argument("--min_lbf", type = "double", default = 2, help = "Minimum LBF value (default: 2)")
  parser$add_argument("--max_p_strict", type = "double", default = 5e-8, help = "Maximum p-value (default: 5e-8)")
  parser$add_argument("--max_p_lenient", type = "double", default = 1e-5, help = "Maximum p-value (default: 1e-5)")
  parser$add_argument("--output", required = TRUE, help = "Name of the output file")

  # Arguments
  args <- parser$parse_args(argv)

  # Toon de waarden (optioneel, voor debug)
  print(args)

  # Gebruik de argumenten in je script
  loci_file <- args$loci
  min_lbf <- args$min_lbf
  max_p_strict <- as.numeric(args$max_p_strict)
  max_p_lenient <- as.numeric(args$max_p_lenient)
  output_file <- args$output

  # Voorbeeld: print de waarden
  cat("Loci bestand:", loci_file, "\n")
  cat("Min LBF:", min_lbf, "\n")
  cat("Max P strict:", max_p_strict, "\n")
  cat("Max P lenient:", max_p_lenient, "\n")
  cat("Output bestand:", output_file, "\n")

  # Function to process a single file
  process_file <- function(file, min_lbf, max_p) {
    dt <- fread(file)

    # Compute p-value
    dt[, p_value := 2 * pnorm(abs(beta/se), lower.tail = FALSE)]

    # Identify all LBF columns
    lbf_cols <- grep("^lbf_cs_", names(dt), value = TRUE)

    # Check each LBF column
    summary_list <- lapply(lbf_cols, function(col) {

      # Only lbf
      has_signal_lbf <- any(dt[[col]] >= min_lbf)
      # Strict
      has_signal_lbf_pstrict <- any(dt[[col]] >= min_lbf & dt$p_value <= max_p_strict)
      # Lenient
      has_signal_lbf_plenient <- any(dt[[col]] >= min_lbf & dt$p_value <= max_p_lenient)

      # List results
      list(LBF_Column = col,
           Has_Signal_Lbf_P_strict = has_signal_lbf_pstrict,
           Has_Signal_Lbf_P_lenient = has_signal_lbf_plenient,
           Has_Signal_Lbf = has_signal_lbf)
    })

    # Combine into a data.table
    summary_dt <- rbindlist(summary_list)
    summary_dt[, File := basename(file)]  # Add filename
    return(summary_dt)
  }

  # Apply to all loci files using lapply and combine results
  all_summary <- rbindlist(lapply(args$loci, process_file, min_lbf = args$min_lbf, max_p = args$max_p))

  # Write output
  fwrite(all_summary, args$output, row.names=F, col.names=T, sep="\t")
  cat("Summary written to", args$output, "\n")

  # Process output
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}