#!/usr/bin/env Rscript

library(coloc)
library(arrow)
library(data.table)
library(stringr)
library(dplyr)
library(argparse)
library(duckdbfs)

setDTthreads(1)

RES_COL_NAMES <- c(
  "ctc_gwas", "ctc_locus", "GWAS", "GWAS_locus", "nsnps", "hit1", "hit2",
  "signal1", "signal2", "PP_H0", "PP_H1", "PP_H2", "PP_H3", "PP_H4", "coloc_cs95_size",
  "cs95_variants", "cs95_PIP", "eQTL_lead_eQTL_Z", "eQTL_lead_GWAS_Z"
)

parser <- ArgumentParser(description = "Run coloc analysis.")

parser$add_argument("--overlap_file",
  metavar = "file",
  type = "character",
  help = "File containing eqtl-gwas region overlaps."
)

parser$add_argument("--coloc_threshold",
  type = "numeric",
  help = "What PP4 threshold to use to declare colocalisation."
)

parser$add_argument("--full_results",
  type = "logical",
  help = "Whether to output all coloc results, including those which show no colocalisation."
)

parser$add_argument("--reference",
  metavar = "file", type = "character",
  help = "eQTLGen SNP reference file in parquet format."
)

parser$add_argument("--min_lbf_threshold",
                    metavar = "float", type = "numeric", default=2.0,
                    help = "Minimum LBF threshold for GWAS credible sets (exclusive)"
)

args <- parser$parse_args()

args$full_results[args$full_results == "true"] <- TRUE
args$full_results[args$full_results == "false"] <- FALSE

args$mrlink2[args$mrlink2 == "true"] <- TRUE
args$mrlink2[args$mrlink2 == "false"] <- FALSE

overlaps <- fread(args$overlap_file)
ctc_traits <- unique(overlaps$gwas_a)
gwas_files <- unique(overlaps$b_file)

min_lbf_threshold <- args$min_lbf_threshold

# Create a summary table. In essence this could be the table in which we write all colocalization results,
# But let's try to use this just for debugging purposses.
# - Is the overlap eligeble for Coloc after QC
# - Is  the overlap show colocalization
results_summary <- data.frame(overlaps, stringsAsFactors = F)
results_summary$pass_qc <- NA
results_summary$n_ctc_lbf_vectors_before_filtering <- NA
results_summary$n_ctc_lbf_vectors_after_filtering <- NA
results_summary$n_gwas_lbf_vectors_before_filtering <- NA
results_summary$n_gwas_lbf_vectors_after_filtering <- NA
results_summary$n_colocalizing_signal_pairs <- NA
results_summary$n_attempted_signal_pairs <- NA
results_summary$n_signal_combinations <- NA

str(results_summary)

fwrite(overlaps, "overlaps.txt")

ctc_gwas_iterator <- 0

message(paste(nrow(overlaps), "COLOC analyses"))

indik <- 1
for (ctc_gwas_indik in ctc_traits) {
  # Make output table
  res <- data.table(
    ctc_gwas  = NA,
    ctc_locus = NA,
    GWAS = NA,
    GWAS_locus = NA,
    nsnps = NA,
    hit1 = NA,
    hit2 = NA,
    signal1 = NA,
    signal2 = NA,
    PP_H0 = NA,
    PP_H1 = NA,
    PP_H2 = NA,
    PP_H3 = NA,
    PP_H4 = NA,
    coloc_cs95_size = NA,
    cs95_variants = NA,
    cs95_PIP = NA,
    ctc_lead_ctc_Z = NA,
    ctc_lead_GWAS_Z = NA
  )[-1]

  fwrite(res, paste0(ctc_gwas_indik, "_coloc_results.txt"), sep = "\t")

  ctc_gwas_iterator <- ctc_gwas_iterator + 1

  inp_coloc <- overlaps[gwas_a == ctc_gwas_indik]
  ctc_files <- unique(inp_coloc$a_file)

  ctc_locus_iterator <- 0
  ctc_locus_length <- length(ctc_files)
  for (ctc_gwas_file in ctc_files) {
    ctc_locus_iterator <- ctc_locus_iterator + 1
    message(paste("eQTL:", ctc_gwas_file))
    inp_coloc_f <- inp_coloc[a_file == ctc_gwas_file]
    message(paste(nrow(inp_coloc_f), "overlapping GWAS loci."))

    eqtl <- fread(file.path("GWAS_finemap", ctc_gwas_file))
    eqtl <- eqtl[order(variant_index)]

    eqtl_variant_indices <- eqtl$variant_index

    eqtl_mat <- as.matrix(eqtl[, str_detect(colnames(eqtl), "lbf_cs_"), with = FALSE])

    # Filter CTC gwas mat on lbf vectors with at least an LBF > 2
    ctc_cs_lbf_pass <- apply(eqtl_mat, 2, max) > min_lbf_threshold

    # Filter on credible sets where the minimal LBF threshold is at least more than 2
    if (!any(ctc_cs_lbf_pass)) {
      print(sprintf("No CTC GWAS credible sets left after selecting LBF vectors with at least an LBF greater than %f", min_lbf_threshold))
      next
    }

    eqtl_mat <- eqtl_mat[ ,ctc_cs_lbf_pass, drop=F]

    eqtl_signals <- as.integer(str_remove(colnames(eqtl_mat)[colSums(eqtl_mat) != 0 & !is.na(colSums(eqtl_mat))], "lbf_cs_"))

    # SNP reference
    ref <- duckdbfs::open_dataset(args$reference) %>%
      filter(variant_index %in% eqtl$variant_index) %>%
      select(variant_index, variant) %>%
      collect() %>%
      as.data.table()

    ref <- ref[order(variant_index), ]

    rownames(eqtl_mat) <- ref$variant

    gwas_loci_iterator <- 0
    length_gwas_loci <- nrow(inp_coloc_f)

    message("Iterating over GWASs...")

    # colocalisation
    for (gwas_ind in seq_len(nrow(inp_coloc_f))) {
      focal_gwas_locus_file <- inp_coloc_f$b_file[gwas_ind]
      message(paste("Analysing: ", focal_gwas_locus_file))

      print(focal_gwas_locus_file)

      overlap_row <- results_summary$a_file == ctc_gwas_file & results_summary$b_file == focal_gwas_locus_file
      results_summary[overlap_row, "n_ctc_lbf_vectors_before_filtering"] <- length(ctc_cs_lbf_pass)
      results_summary[overlap_row, "n_ctc_lbf_vectors_after_filtering"] <- sum(ctc_cs_lbf_pass)

      gwas_loci_iterator <- gwas_loci_iterator + 1
      gwas <- fread(file.path("GWAS_finemap", focal_gwas_locus_file))

      gwas <- gwas[order(variant_index)]
      gwas <- gwas[variant_index %in% eqtl$variant_index]

      ref2 <- ref[variant_index %in% gwas$variant_index, ]
      gwas_variant_indices <- gwas$variant_index

      gwas_mat <- as.matrix(gwas[, str_detect(colnames(gwas), "lbf_cs_"), with = FALSE])

      # Filter GWAS mat on lbf vectors with at least an LBF > 2
      gwas_cs_lbf_pass <- apply(gwas_mat, 2, max) > min_lbf_threshold

      overlap_row <- results_summary$a_file == ctc_gwas_file & results_summary$b_file == focal_gwas_locus_file
      results_summary[overlap_row, "n_gwas_lbf_vectors_before_filtering"] <- length(gwas_cs_lbf_pass)
      results_summary[overlap_row, "n_gwas_lbf_vectors_after_filtering"] <- sum(gwas_cs_lbf_pass)

      # Filter on credible sets where the minimal LBF threshold is at least more than 2
      if (!any(gwas_cs_lbf_pass)) {
        print(sprintf("No GWAS credible sets left after selecting LBF vectors with at least an LBF greater than %f", min_lbf_threshold))
        next
      }

      gwas_mat <- gwas_mat[ ,gwas_cs_lbf_pass, drop=F]

      rownames(gwas_mat) <- ref2$variant
      rm(ref2)

      gwas_signals <- as.integer(str_remove(colnames(gwas_mat)[colSums(gwas_mat) != 0 & !is.na(colSums(gwas_mat))], "lbf_cs_"))

      eqtl2 <- eqtl[variant_index %in% gwas_variant_indices, ]
      gwas2 <- gwas[variant_index %in% eqtl_variant_indices, ]

      passes_threshold <- nrow(gwas2) > 350 & length(gwas_signals) > 0
      results_summary[overlap_row, "n_signal_combinations"] <- length(eqtl_signals) * length(gwas_signals)
      results_summary[overlap_row, "pass_qc"] <- passes_threshold

      n_colocalizating_signal_pairs <- 0
      n_attempted_signal_pairs <- 0

      if (passes_threshold) {
        message("Trying analysis!")
        message("Iterating over independent eQTL signals...")
        for (i in eqtl_signals) {
          message("Iterating over independent GWAS signals...")
          for (j in gwas_signals) {
            print(sprintf("eQTL signal %s, GWAS signal %s", i, j))
            res_coloc <- coloc.bf_bf(
              eqtl_mat[, paste0("lbf_cs_", i)],
              gwas_mat[, paste0("lbf_cs_", j)],
              p1 = 1e-04,
              p2 = 1e-04,
              overlap.min = 0.5,
              trim_by_posterior = TRUE
            )

            if ("results" %in% names(res_coloc)) {
              n_attempted_signal_pairs <- n_attempted_signal_pairs + 1
            }

            gc()

            res_temp <- data.table(
              ctc_gwas = inp_coloc_f$gwas_a[gwas_ind],
              ctc_locus = paste0(inp_coloc_f$loc_chr[gwas_ind], ":", inp_coloc_f$a_loc_start[gwas_ind], "-", inp_coloc_f$a_loc_end[gwas_ind]),
              GWAS = inp_coloc_f$gwas_b[gwas_ind],
              GWAS_locus = paste0(inp_coloc_f$loc_chr[gwas_ind], ":", inp_coloc_f$b_loc_start[gwas_ind], "-", inp_coloc_f$b_loc_end[gwas_ind]),
              res_coloc$summary
            )

            res_temp$idx1 <- i
            res_temp$idx2 <- j

            pass_coloc <- !is.null(res_temp$PP.H4.abf) &&
              !is.na(res_temp$PP.H4.abf) &&
              length(res_temp$PP.H4.abf) > 0 &&
              res_temp$PP.H4.abf > args$coloc_threshold

            n_colocalizating_signal_pairs <- n_colocalizating_signal_pairs + as.integer(pass_coloc)

            if (gwas_ind == 1 & i == min(eqtl_signals) & j == min(gwas_signals)) {
              if (pass_coloc) {
                o <- order(res_coloc$results$SNP.PP.H4, decreasing = TRUE)
                cs <- cumsum(res_coloc$results$SNP.PP.H4[o])
                w <- which(cs > 0.95)[1]
                res_temp$coloc_cs_size <- length(res_coloc$results[o, ][1:w, ]$snp)
                res_temp$cs_variants <- paste(res_coloc$results[o, ][1:w, ]$snp, collapse = "; ")
                res_temp$cs_PP <- paste(res_coloc$results[o, ][1:w, ]$`SNP.PP.H4.abf`, collapse = "; ")
                res_temp$ctc_Z <- eqtl2[variant_index == ref[variant == res_temp$hit1]$variant_index]$beta / eqtl2[variant_index == ref[variant == res_temp$hit1]$variant_index]$se
                res_temp$GWAS_Z <- gwas[variant_index == ref[variant == res_temp$hit1]$variant_index]$beta / gwas[variant_index == ref[variant == res_temp$hit1]$variant_index]$se
              } else {
                res_temp$coloc_cs_size <- NA
                res_temp$cs_variants <- NA
                res_temp$cs_PP <- NA
                res_temp$ctc_Z <- NA
                res_temp$GWAS_Z <- NA
              }
            } else {
              if (pass_coloc) {
                o <- order(res_coloc$results$SNP.PP.H4, decreasing = TRUE)
                cs <- cumsum(res_coloc$results$SNP.PP.H4[o])
                w <- which(cs > 0.95)[1]
                res_temp$coloc_cs_size <- length(res_coloc$results[o, ][1:w, ]$snp)
                res_temp$cs_variants <- paste(res_coloc$results[o, ][1:w, ]$snp, collapse = ";")
                res_temp$cs_PP <- paste(res_coloc$results[o, ][1:w, ]$`SNP.PP.H4.abf`, collapse = ";")
                res_temp$ctc_Z <- eqtl2[variant_index == ref[variant == res_temp$hit1]$variant_index]$beta / eqtl2[variant_index == ref[variant == res_temp$hit1]$variant_index]$se
                res_temp$GWAS_Z <- gwas[variant_index == ref[variant == res_temp$hit1]$variant_index]$beta / gwas[variant_index == ref[variant == res_temp$hit1]$variant_index]$se
              } else {
                res_temp$coloc_cs_size <- NA
                res_temp$cs_variants <- NA
                res_temp$cs_PP <- NA
                res_temp$ctc_Z <- NA
                res_temp$GWAS_Z <- NA
              }
            }
            if (ncol(res_temp) == length(RES_COL_NAMES)) {
              res_temp <- res_temp[, c(1:7, 13:14, 8:12, 15:19), with = FALSE]

              colnames(res_temp) <- RES_COL_NAMES
              if (isFALSE(args$full_results)) {
                res_temp <- res_temp[PP_H4 > args$coloc_threshold]
              }

              if (nrow(res_temp) > 0) {
                fwrite(res_temp, paste0(ctc_gwas_indik, "_coloc_results.txt"), sep = "\t", append = TRUE)
              message("Output written!")
              }
              rm(res_temp)
              gc()
              indik <- indik + 1
            }
          }
          message("Iterating over independent GWAS signals...done!")
        }
        message("Iterating over independent eQTL signals...done!")
        message(paste0(
          "Analysed: ", ctc_gwas_iterator, "/", length(ctc_files), " eQTL gene, ",
          ctc_locus_iterator, "/", ctc_locus_length, " eQTL locus, ",
          gwas_loci_iterator, "/", length_gwas_loci, " overlapping GWAS locus"
        ))
        results_summary[overlap_row, "n_colocalizing_signal_pairs"] <- n_colocalizating_signal_pairs
        results_summary[overlap_row, "n_attempted_signal_pairs"] <- n_attempted_signal_pairs

      }
      rm(gwas, gwas2, eqtl2)
      gc()
    }
    rm(eqtl, eqtl_mat)
    gc()
    message(paste("Finalised with eQTL iterations"))
  }
}

fwrite(results_summary, "results_summary_per_overlap.txt", sep = "\t", row.names=F, col.names=T, quote=F, na="NA")

warnings()
