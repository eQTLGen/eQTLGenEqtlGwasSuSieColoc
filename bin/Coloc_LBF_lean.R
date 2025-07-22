#!/usr/bin/env Rscript

library(coloc)
library(arrow)
library(data.table)
library(stringr)
library(dplyr)
library(argparse)
library(duckdbfs)

setDTthreads(1)
parser <- ArgumentParser(description = "Run coloc analysis.")

parser$add_argument("--overlap_file",
  metavar = "file",
  type = "character",
  help = "File containing eqtl-gwas region overlaps."
)

parser$add_argument("--gene_names",
  metavar = "file",
  type = "character",
  help = "File containing ENSG and HGNC SYMBOL mappings."
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

args <- parser$parse_args()

args$full_results[args$full_results == "true"] <- TRUE
args$full_results[args$full_results == "false"] <- FALSE

args$mrlink2[args$mrlink2 == "true"] <- TRUE
args$mrlink2[args$mrlink2 == "false"] <- FALSE

overlaps <- fread(args$overlap_file)

eqtl_genes <- list.files()[str_detect(list.files(), "cis") | str_detect(list.files(), "trans")]
eqtl_genes <- unique(str_remove(eqtl_genes, "__.*"))

message("eQTL genes:")
print(eqtl_genes)

gwas_files <- list.files()[str_detect(list.files(), "gwas")]

overlaps <- overlaps[eqtl_gene %in% eqtl_genes & gwas_file %in% gwas_files]

fwrite(overlaps, "overlaps.txt")

eqtl_gene_iterator <- 0

message(paste(nrow(overlaps), "COLOC analyses"))

indik <- 1
for (eqtl_gene_indik in eqtl_genes) {
  # Make output table
  res <- data.table(
    eQTL_gene = NA,
    eQTL_locus = NA,
    eQTL_type = NA,
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
    eQTL_lead_eQTL_Z = NA,
    eQTL_lead_GWAS_Z = NA
  )[-1]

  fwrite(res, paste0(eqtl_gene_indik, "_coloc_results.txt"), sep = "\t")

  eqtl_gene_iterator <- eqtl_gene_iterator + 1

  inp_coloc <- overlaps[eqtl_gene == eqtl_gene_indik]
  eqtl_file <- unique(inp_coloc$eqtl_file)

  eqtl_locus_iterator <- 0
  eqtl_locus_length <- length(eqtl_file)
  for (eqtl_locus in eqtl_file) {
    eqtl_locus_iterator <- eqtl_locus_iterator + 1
    message(paste("eQTL:", eqtl_locus))
    inp_coloc_f <- inp_coloc[eqtl_file == eqtl_locus]
    message(paste(nrow(inp_coloc_f), "overlapping GWAS loci."))

    eqtl <- fread(eqtl_locus)
    eqtl <- eqtl[order(variant_index)]

    eqtl_variant_indices <- eqtl$variant_index

    eqtl_mat <- as.matrix(eqtl[, str_detect(colnames(eqtl), "lbf_cs_"), with = FALSE])

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
    for (gwas_ind in 1:nrow(inp_coloc_f)) {
      message(paste("Analysing: ", inp_coloc_f$gwas_file[gwas_ind]))

      print(inp_coloc_f$gwas_file[gwas_ind])

      gwas_loci_iterator <- gwas_loci_iterator + 1
      gwas <- fread(inp_coloc_f$gwas_file[gwas_ind])

      gwas <- gwas[order(variant_index)]
      gwas <- gwas[variant_index %in% eqtl$variant_index]

      ref2 <- ref[variant_index %in% gwas$variant_index, ]
      gwas_variant_indices <- gwas$variant_index

      gwas_mat <- as.matrix(gwas[, str_detect(colnames(gwas), "lbf_cs_"), with = FALSE])

      rownames(gwas_mat) <- ref2$variant
      rm(ref2)

      gwas_signals <- as.integer(str_remove(colnames(gwas_mat)[colSums(gwas_mat) != 0 & !is.na(colSums(gwas_mat))], "lbf_cs_"))

      eqtl2 <- eqtl[variant_index %in% gwas_variant_indices, ]
      gwas2 <- gwas[variant_index %in% eqtl_variant_indices, ]

      if (nrow(gwas2) > 500 & length(gwas_signals) > 0) {
        message("Trying analysis!")
        message("Iterating over independent eQTL signals...")
        for (i in eqtl_signals) {
          message("Iterating over independent GWAS signals...")
          for (j in gwas_signals) {
            res_coloc <- coloc.bf_bf(
              eqtl_mat[, i],
              gwas_mat[, j],
              p1 = 1e-04,
              p2 = 1e-04,
              overlap.min = 0.5,
              trim_by_posterior = TRUE
            )

            gc()

            res_temp <- data.table(
              eQTL_gene = inp_coloc_f$eqtl_gene[gwas_ind],
              eQTL_locus = paste0(inp_coloc_f$loc_chr[gwas_ind], ":", inp_coloc_f$eqtl_loc_start[gwas_ind], "-", inp_coloc_f$eqtl_loc_end[gwas_ind]),
              eQTL_type = str_remove(inp_coloc_f$eqtl_file[gwas_ind], ".*___"),
              GWAS = inp_coloc_f$gwas[gwas_ind],
              GWAS_locus = paste0(inp_coloc_f$loc_chr[gwas_ind], ":", inp_coloc_f$gwas_loc_start[gwas_ind], "-", inp_coloc_f$gwas_loc_end[gwas_ind]),
              res_coloc$summary
            )

            res_temp$eQTL_type <- str_remove(res_temp$eQTL_type, ".txt.gz")

            res_temp$idx1 <- i
            res_temp$idx2 <- j


            if (gwas_ind == 1 & i == min(eqtl_signals) & j == min(gwas_signals)) {
              if (!is.null(res_temp$PP.H4.abf) &&
                !is.na(res_temp$PP.H4.abf) &&
                length(res_temp$PP.H4.abf) > 0 &&
                res_temp$PP.H4.abf > args$coloc_threshold) {
                o <- order(res_coloc$results$SNP.PP.H4, decreasing = TRUE)
                cs <- cumsum(res_coloc$results$SNP.PP.H4[o])
                w <- which(cs > 0.95)[1]
                res_temp$coloc_cs_size <- length(res_coloc$results[o, ][1:w, ]$snp)
                res_temp$cs_variants <- paste(res_coloc$results[o, ][1:w, ]$snp, collapse = "; ")
                res_temp$cs_PP <- paste(res_coloc$results[o, ][1:w, ]$`SNP.PP.H4.abf`, collapse = "; ")
                res_temp$eQTL_Z <- eqtl2[variant_index == ref[variant == res_temp$hit1]$variant_index]$beta / eqtl2[variant_index == ref[variant == res_temp$hit1]$variant_index]$standard_error
                res_temp$GWAS_Z <- gwas[variant_index == ref[variant == res_temp$hit1]$variant_index]$beta / gwas[variant_index == ref[variant == res_temp$hit1]$variant_index]$se
              } else {
                res_temp$coloc_cs_size <- NA
                res_temp$cs_variants <- NA
                res_temp$cs_PP <- NA
                res_temp$eQTL_Z <- NA
                res_temp$GWAS_Z <- NA
              }
            } else {
              if (!is.null(res_temp$PP.H4.abf) &&
                !is.na(res_temp$PP.H4.abf) &&
                length(res_temp$PP.H4.abf) > 0 &&
                res_temp$PP.H4.abf > args$coloc_threshold) {
                o <- order(res_coloc$results$SNP.PP.H4, decreasing = TRUE)
                cs <- cumsum(res_coloc$results$SNP.PP.H4[o])
                w <- which(cs > 0.95)[1]
                res_temp$coloc_cs_size <- length(res_coloc$results[o, ][1:w, ]$snp)
                res_temp$cs_variants <- paste(res_coloc$results[o, ][1:w, ]$snp, collapse = ";")
                res_temp$cs_PP <- paste(res_coloc$results[o, ][1:w, ]$`SNP.PP.H4.abf`, collapse = ";")
                res_temp$eQTL_Z <- eqtl2[variant_index == ref[variant == res_temp$hit1]$variant_index]$beta / eqtl2[variant_index == ref[variant == res_temp$hit1]$variant_index]$standard_error
                res_temp$GWAS_Z <- gwas[variant_index == ref[variant == res_temp$hit1]$variant_index]$beta / gwas[variant_index == ref[variant == res_temp$hit1]$variant_index]$se
              } else {
                res_temp$coloc_cs_size <- NA
                res_temp$cs_variants <- NA
                res_temp$cs_PP <- NA
                res_temp$eQTL_Z <- NA
                res_temp$GWAS_Z <- NA
              }
            }
            if (ncol(res_temp) == 20) {
              res_temp <- res_temp[, c(1:8, 14:15, 9:13, 16:20), with = FALSE]
              colnames(res_temp) <- c(
                "eQTL_gene", "eQTL_locus", "eQTL_type", "GWAS", "GWAS_locus", "nsnps", "hit1", "hit2",
                "signal1", "signal2", "PP_H0", "PP_H1", "PP_H2", "PP_H3", "PP_H4", "coloc_cs95_size",
                "cs95_variants", "cs95_PIP", "eQTL_lead_eQTL_Z", "eQTL_lead_GWAS_Z"
              )
              if (isFALSE(args$full_results)) {
                res_temp <- res_temp[PP_H4 > args$coloc_threshold]
              }

              if (nrow(res_temp) > 0) {
                fwrite(res_temp, paste0(eqtl_gene_indik, "_coloc_results.txt"), sep = "\t", append = TRUE)
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
          "Analysed: ", eqtl_gene_iterator, "/10 eQTL gene, ",
          eqtl_locus_iterator, "/", eqtl_locus_length, " eQTL locus, ",
          gwas_loci_iterator, "/", length_gwas_loci, " overlapping GWAS locus"
        ))
      }
      rm(gwas, gwas2, eqtl2)
      gc()
    }
    rm(eqtl, eqtl_mat)
    gc()
    message(paste("Finalised with eQTL iterations"))
  }
}
warnings()
