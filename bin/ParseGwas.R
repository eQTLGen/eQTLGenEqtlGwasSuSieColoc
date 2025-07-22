#!/usr/bin/env Rscript

library(coloc)
library(arrow)
library(data.table)
library(stringr)
library(dplyr)
library(IGUtilityPackage)
library(argparse)

set.seed(8)

setDTthreads(1)
parser <- ArgumentParser(description = "Parse GWAS summary statistics.")

parser$add_argument("--gwas_folder",
  metavar = "file",
  type = "character",
  help = "Folder containing GWAS parquet files."
)

parser$add_argument("--gwas_id",
  type = "character",
  help = "GWAS ID for which to run the analysis."
)

parser$add_argument("--win_size",
  type = "numeric",
  help = "GWAS window size for defining loci."
)

args <- parser$parse_args()

# functions
FindLoci <- function(dt, p_thresh = 5e-8, window_size = 500000) {
  # Ensure numeric p-values
  dt[, p := as.numeric(p)]

  # Identify significant variants
  sig <- dt[p < p_thresh]

  if (nrow(sig) == 0) {
    return(data.table()) # No significant loci
  }

  # Define ±window windows
  sig_windows <- sig[, .(chromosome, start = bp - window_size, end = bp + window_size)]

  # Sort windows for merging
  setorder(sig_windows, chromosome, start)
  merged <- sig_windows[1]

  # Progress bar setup
  pb <- txtProgressBar(min = 1, max = nrow(sig_windows), style = 3)

  for (i in 2:nrow(sig_windows)) {
    setTxtProgressBar(pb, i)

    cur_chr <- sig_windows[i, chromosome]
    cur_start <- sig_windows[i, start]
    cur_end <- sig_windows[i, end]

    if (cur_chr == merged[.N, chromosome] && cur_start <= merged[.N, end]) {
      merged[.N, end := max(merged[.N, end], cur_end)]
    } else {
      merged <- rbind(merged, sig_windows[i])
    }
  }
  close(pb)

  # Create locus_id
  merged[, locus_id := paste0(chromosome, ":", start, "_", end)]

  # Count variants per locus
  setkey(dt, chromosome, bp)
  loci_summary <- merged[,
    {
      nvar <- dt[chromosome == .BY$chromosome & bp >= start & bp <= end, .N]
      .(start, end, locus_id, n_variants = nvar)
    },
    by = .(chromosome)
  ]

  return(loci_summary)
}

IdentifyLeadSNPs <- function(data,
                             window = 1000000,
                             Pthresh = 5e-8,
                             snp_id_col = "snp",
                             snp_chr_col = "chr",
                             snp_pos_col = "pos",
                             eff_all_col = "ea",
                             other_all_col = "nea",
                             beta_col = "beta",
                             se_col = "se",
                             p_col = NULL,
                             loci = FALSE,
                             merge_overlapping = FALSE,
                             verbose = TRUE) {
  if (is.null(p_col) & !is.null(beta_col) & !is.null(se_col)) {
    data <- data.table(
      SNP = data[[snp_id_col]],
      chr = data[[snp_chr_col]],
      pos = data[[snp_pos_col]],
      ea = data[[eff_all_col]],
      nea = data[[other_all_col]],
      beta = data[[beta_col]],
      se = data[[se_col]]
    )

    data$P <- ZtoP(data$beta / data$se)
    data$sig_ind <- data$beta / data$se
  } else if (!is.null(p_col) & is.null(beta_col) & is.null(se_col)) {
    data <- data.table(
      SNP = data[[snp_id_col]],
      chr = data[[snp_chr_col]],
      pos = data[[snp_pos_col]],
      ea = data[[eff_all_col]],
      nea = data[[other_all_col]],
      P = data[[p_col]]
    )

    data$sig_ind <- -log10(P)
  } else if (!is.null(p_col) & !is.null(beta_col) & !is.null(se_col)) {
    data <- data.table(
      SNP = data[[snp_id_col]],
      chr = data[[snp_chr_col]],
      pos = data[[snp_pos_col]],
      ea = data[[eff_all_col]],
      nea = data[[other_all_col]],
      beta = data[[beta_col]],
      se = data[[se_col]],
      P = data[[p_col]]
    )
    data$sig_ind <- data$beta / data$se
  } else if (is.null(p_col) & is.null(beta_col) & is.null(se_col)) {
    message("There is no beta/se nor P-value in the data")
    stop
  }

  data_f <- data[data$P < as.numeric(Pthresh), ]

  # Iteratively identify most significant SNP, and remove all other SNPs in the window
  res <- data_f[-c(1:nrow(data_f)), ]

  while (min(data_f$P) <= Pthresh) {
    lead_snp <- data_f[abs(data_f$sig_ind) == max(abs(data_f$sig_ind)), ]
    if (nrow(lead_snp) > 1) {
      lead_snp <- lead_snp[1, ]
    }
    res <- rbind(res, lead_snp)
    data_f <- data_f[!(data_f$chr == lead_snp$chr & data_f$pos > lead_snp$pos - window & data_f$pos < lead_snp$pos + window), ]
    if (isTRUE(verbose)) {
      message(paste0("Added: ", lead_snp$SNP, " | ", lead_snp$chr, ":", lead_snp$pos))
    }
    if (nrow(data_f) == 0) {
      break
    }
  }

  res <- res[order(chr, pos), -ncol(res), with = FALSE]

if (isTRUE(loci)) {
  message("Extracting full loci...")

  temp_list <- vector("list", nrow(res))

  for (i in seq_len(nrow(res))) {
    chr_i <- res[i, chr]
    pos_i <- res[i, pos]
    snp_i <- res[i, SNP]

    # Extract window
    locus <- data[chr == chr_i & pos > (pos_i - window) & pos < (pos_i + window), .SD, .SDcols = !ncol(data)]
    locus[, `:=`(locus = i, lead = SNP == snp_i)]
    temp_list[[i]] <- locus

    if (isTRUE(verbose)) {
      message(paste0(i, "/", nrow(res)))
    }
  }

  temp_res <- rbindlist(temp_list)

  if (isTRUE(merge_overlapping)) {
    message("Merging overlapping loci...")

    # Create windows table
    win_dt <- res[, .(chr, start = pos - window, end = pos + window)]

    # Sort and merge overlapping windows
    setorder(win_dt, chr, start)
    merged <- win_dt[1]

    for (i in 2:nrow(win_dt)) {
      if (win_dt[i, chr] == merged[.N, chr] && win_dt[i, start] <= merged[.N, end]) {
        merged[.N, end := max(merged[.N, end], win_dt[i, end])]
      } else {
        merged <- rbind(merged, win_dt[i])
      }
      message(paste0(i, "/", nrow(win_dt)))
    }

    # Assign merged locus IDs
    temp_res[, merged_locus := NA_integer_]
    for (i in seq_len(nrow(merged))) {
      temp_res[chr == merged[i, chr] & pos >= merged[i, start] & pos <= merged[i, end], merged_locus := i]
    }

    # Replace original locus column with merged_locus if merging is enabled
    temp_res[, locus := merged_locus][, merged_locus := NULL]
  }

  res <- temp_res

}
  return(res)
}




print(args)

gwas <- open_dataset(paste0(args$gwas_folder, "/gwas_id=", args$gwas_id)) %>%
  collect() %>%
  as.data.table()

message(paste("Rows in data:", nrow(gwas)))

message("Removing rows with beta=0 and/or se=0...")

gwas <- gwas[beta != 0 & se != 0]

message("Removing rows with beta=0 and/or se=0...done!")

message("Removing HLA region...")
gwas <- gwas[!(chromosome == 6 & bp > 28510120 & bp < 33480577)]
message("Removing HLA region...done!")

message(paste("Rows in data:", nrow(gwas)))

message("Finding loci...")
Loci <- IdentifyLeadSNPs(gwas, 
snp_id_col = "variant_index", 
snp_chr_col = "chromosome", 
snp_pos_col = "bp", 
eff_all_col = "str_allele2", 
other_all_col = "str_allele1", 
beta_col = "beta", 
se_col = "se", 
loci = FALSE,
window = 1000000)

Loci <- Loci %>% 
rowwise() %>% 
summarise(chromosome = unique(chr), start = pos - 1000000, end = pos + 1000000, locus_id = paste0(chr, ":", start, "_", end)) %>% 
as.data.table()
message("Finding loci...done!")

rm(gwas)
gc()

gwas <- open_dataset(paste0(args$gwas_folder, "/gwas_id=", args$gwas_id))

nr_loci <- nrow(Loci)
message(paste(nr_loci, "loci found!"))

for (i in 1:nrow(Loci)) {

  message(paste0("Analysing: ", Loci$locus_id[i], " ", i, "/", nrow(Loci), "..."))

  start.time <- Sys.time()

  temp_chr <- as.integer(Loci$chromosome[i])
  temp_start <- as.integer(Loci$start[i])
  temp_end <- as.integer(Loci$end[i])

  locus <- gwas %>%
    filter(chromosome == temp_chr & bp > temp_start & bp < temp_end) %>%
    select(variant_index, chromosome, bp, beta, se, p, N) %>%
    collect() %>%
    filter(N >= 0.8 * max(N)) %>%
    as.data.table()

  message(paste(nrow(locus), "variants"))

  fwrite(locus, paste0(args$gwas_id, "__", Loci$locus_id[i], "___", "temp", ".txt.gz"), sep = "\t")

  end.time <- Sys.time()
  message(paste("Elapsed:", end.time - start.time))
  message(paste0("Analysing: ", Loci$locus_id[i], " ", i, "/", nrow(Loci), "...done!"))
  gc()
}

warnings()
