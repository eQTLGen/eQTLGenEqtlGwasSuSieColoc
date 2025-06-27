#!/usr/bin/env Rscript

library(coloc)
library(arrow)
library(data.table)
library(stringr)
library(dplyr)
library(argparse)

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

parser$add_argument("--mrlink2",
  type = "logical",
  help = "Whether to run also MRLink2."
)

parser$add_argument("--reference",
  metavar = "file", type = "character",
  help = "eQTLGen SNP reference file in parquet format."
)

parser$add_argument("--ld_folder",
  metavar = "file",
  help = "Folder containing permuted LD files."
)

args <- parser$parse_args()

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

# MRLink2 functions from Dr. Adriaan van der Graaf (https://github.com/adriaan-vd-graaf/mrlink2/blob/main/R/mr_link_2_functions.R)
mr_link2_loglik_reference_v2 <- function(th, lam, c_x, c_y, n_x, n_y) {
  # Convert n_x and n_y to float
  n_x <- as.numeric(n_x)
  n_y <- as.numeric(n_y)
  a <- th[1]
  tX <- abs(th[2])
  tY <- abs(th[3])

  Dyy <- 1 / (n_y * lam + tY)

  if (a != 0) {
    Dxx <- 1 / (exp(log(a^2 * n_y + n_x) + log(lam)) + tX -
      exp(log(a^2 * n_y^2 * lam^2) - log(n_y * lam + tY)))
    Dxy <- -Dxx * a * exp(log(n_y * lam) - log(n_y * lam + tY))
    Dyy <- Dyy + exp(log(Dxx * (a^2 * n_y^2 * lam^2)) - (2 * log(n_y * lam + tY)))
    asq_ny_sq_lam_sq_div_ny_lam_ty <- exp(log(a^2 * n_y^2 * lam^2) - log(n_y * lam + tY))
  } else {
    Dxx <- 1 / (exp(log(n_x) + log(lam)) + tX)
    Dxy <- -Dxx * a * exp(log(n_y * lam) - log(n_y * lam + tY))
    asq_ny_sq_lam_sq_div_ny_lam_ty <- 0 * lam
  }

  dX <- n_x * c_x + a * n_y * c_y
  dY <- n_y * c_y
  m <- length(c_x)

  loglik <- -m * log(2 * pi) -
    (1 / 2) * sum(log((a^2 * n_y + n_x) * lam + tX - asq_ny_sq_lam_sq_div_ny_lam_ty)) -
    (1 / 2) * sum(log(n_y * lam + tY)) +
    (1 / 2) * (sum(dX^2 * Dxx) + 2 * sum(dX * dY * Dxy) + sum(dY^2 * Dyy)) -
    (n_x / 2) * sum((c_x^2) / lam) -
    (n_y / 2) * sum((c_y^2) / lam) +
    (m / 2) * (log(n_x) + log(n_y)) - sum(log(lam)) + (m / 2) * (log(tX) + log(tY))

  return(-loglik)
}

mr_link2_loglik_alpha_h0 <- function(th, lam, cX, cY, nX, nY) {
  return(mr_link2_loglik_reference_v2(c(0, th[1], th[2]), lam, cX, cY, nX, nY))
}

mr_link2_loglik_sigma_y_h0 <- function(th, lam, c_x, c_y, n_x, n_y) {
  n_x <- as.numeric(n_x)
  n_y <- as.numeric(n_y)
  a <- th[1]
  tX <- abs(th[2])

  Dyy <- rep(0, length(lam))

  if (a != 0) {
    Dxx <- 1 / (exp(log(a^2 * n_y + n_x) + log(lam)) + tX)
    Dxy <- rep(0, length(lam))
    asq_ny_sq_lam_sq_div_ny_lam_ty <- rep(0, length(lam))
  } else {
    Dxx <- 1 / (exp(log(n_x) + log(lam)) + tX)
    Dxy <- rep(0, length(lam))
    asq_ny_sq_lam_sq_div_ny_lam_ty <- rep(0, length(lam))
  }

  dX <- n_x * c_x + a * n_y * c_y
  dY <- n_y * c_y
  m <- length(c_x)

  loglik <- -m * log(2 * pi) -
    (1 / 2) * sum(log((a^2 * n_y + n_x) * lam + tX - asq_ny_sq_lam_sq_div_ny_lam_ty)) +
    (1 / 2) * (sum(dX^2 * Dxx) + 2 * sum(dX * dY * Dxy) + sum(dY^2 * Dyy)) -
    (n_x / 2) * sum((c_x^2) / lam) -
    (n_y / 2) * sum((c_y^2) / lam) +
    (m / 2) * (log(n_x) + log(n_y)) - sum(log(lam)) + (m / 2) * log(tX)

  return(-loglik)
}



mr_link2 <- function(selected_eigenvalues, selected_eigenvectors,
                     exposure_betas, outcome_betas,
                     n_exp, n_out, sigma_exp_guess, sigma_out_guess) {
  start_time <- Sys.time()

  # Define optimization options
  method <- "Nelder-Mead"
  control <- list(maxit = 300)

  # Calculate c_x and c_y
  c_x <- t(selected_eigenvectors) %*% exposure_betas
  c_y <- t(selected_eigenvectors) %*% outcome_betas

  max_sigma <- sqrt(.Machine$double.xmax)

  # Alpha h0 estimation
  alpha_h0_guesses <- list(
    c(sigma_exp_guess, sigma_out_guess),
    c(max_sigma, max_sigma),
    c(1, max_sigma),
    c(1e3, 1e3)
  )

  alpha_h0_results <- optim(
    par = alpha_h0_guesses[[1]],
    fn = mr_link2_loglik_alpha_h0,
    gr = NULL,
    selected_eigenvalues, c_x, c_y, n_exp, n_out,
    method = method, control = control
  )

  # Loop over other guesses
  for (alpha_h0_guess in alpha_h0_guesses[-1]) {
    if (alpha_h0_results$convergence == 0) break

    new_alpha_h0_results <- optim(
      par = alpha_h0_guess,
      fn = mr_link2_loglik_alpha_h0,
      gr = NULL,
      selected_eigenvalues, c_x, c_y, n_exp, n_out,
      method = method, control = control
    )

    if (alpha_h0_results$value >= new_alpha_h0_results$value) {
      alpha_h0_results <- new_alpha_h0_results
    }
  }

  # Sigma_y estimation
  sigma_y_guesses <- list(
    c(0.0, sigma_exp_guess),
    c(1.0, sigma_exp_guess),
    c(0.0, alpha_h0_results$par[1]),
    c(1.0, alpha_h0_results$par[1]),
    c(0.0, max_sigma),
    c(1e-10, max_sigma)
  )

  sigma_y_h0_results <- optim(
    par = sigma_y_guesses[[1]],
    fn = mr_link2_loglik_sigma_y_h0,
    gr = NULL,
    selected_eigenvalues, c_x, c_y, n_exp, n_out,
    method = method, control = control
  )

  for (sigma_y_guess in sigma_y_guesses[-1]) {
    if (sigma_y_h0_results$convergence == 0) break

    new_sigma_y_h0_results <- optim(
      par = sigma_y_guess,
      fn = mr_link2_loglik_sigma_y_h0,
      gr = NULL,
      selected_eigenvalues, c_x, c_y, n_exp, n_out,
      method = method, control = control
    )

    if (new_sigma_y_h0_results$value < sigma_y_h0_results$value) {
      sigma_y_h0_results <- new_sigma_y_h0_results
    }
  }

  # Ha estimation
  ha_guesses <- list(
    c(0.0, alpha_h0_results$par[1], alpha_h0_results$par[2]),
    c(sigma_y_h0_results$par[1], sigma_y_h0_results$par[2], sqrt(.Machine$double.xmax)),
    c(1.0, alpha_h0_results$par[1], alpha_h0_results$par[2]),
    c(1e-10, max_sigma, max_sigma)
  )

  ha_results <- optim(
    par = ha_guesses[[1]],
    fn = mr_link2_loglik_reference_v2,
    gr = NULL,
    selected_eigenvalues, c_x, c_y, n_exp, n_out,
    method = method, control = control
  )

  for (ha_guess in ha_guesses[-1]) {
    if (ha_results$convergence == 0) break

    new_ha_result <- optim(
      par = ha_guess,
      fn = mr_link2_loglik_reference_v2,
      gr = NULL,
      selected_eigenvalues, c_x, c_y, n_exp, n_out,
      method = method, control = control
    )

    if (new_ha_result$value < ha_results$value) {
      ha_results <- new_ha_result
    }
  }

  # Likelihood Ratio Test and Estimation
  alpha <- ha_results$par[1]
  alpha_chi_sq <- 2 * (alpha_h0_results$value - ha_results$value)
  alpha_p_val <- pchisq(alpha_chi_sq, df = 1, lower.tail = FALSE)
  z_alpha <- ifelse(alpha_chi_sq <= 0, 0.0, sign(alpha) * sqrt(alpha_chi_sq))
  se_alpha <- ifelse(z_alpha != 0, alpha / z_alpha, NA)

  sigma_y <- 1 / abs(ha_results$par[3])
  sigma_y_chi_sq <- 2 * (sigma_y_h0_results$value - ha_results$value)
  sigma_y_p_val <- pchisq(sigma_y_chi_sq, df = 1, lower.tail = FALSE)
  z_sigma_y <- ifelse(sigma_y_chi_sq <= 0, 0.0, sqrt(sigma_y_chi_sq))
  se_sigma_y <- ifelse(z_sigma_y != 0, sigma_y / z_sigma_y, NA)

  # Return results as a list (equivalent to Python's dictionary)
  list(
    alpha = alpha,
    `se(alpha)` = se_alpha,
    `p(alpha)` = alpha_p_val,
    sigma_y = sigma_y,
    `se(sigma_y)` = se_sigma_y,
    `p(sigma_y)` = sigma_y_p_val,
    sigma_x = 1 / abs(ha_results$par[2]),
    alpha_h0_sigma_x = 1 / abs(alpha_h0_results$par[1]),
    alpha_h0_sigma_y = 1 / abs(alpha_h0_results$par[2]),
    alpha_h0_loglik = alpha_h0_results$value,
    sigma_y_h0_alpha = sigma_y_h0_results$par[1],
    sigma_y_h0_sigma_x = 1 / abs(sigma_y_h0_results$par[2]),
    sigma_y_h0_loglik = sigma_y_h0_results$value,
    ha_loglik = ha_results$value,
    optim_alpha_h0_success = alpha_h0_results$convergence == 0,
    optim_alpha_h0_nit = alpha_h0_results$counts,
    optim_sigma_y_h0_success = sigma_y_h0_results$convergence == 0,
    optim_sigma_y_h0_nit = sigma_y_h0_results$counts,
    optim_ha_success = ha_results$convergence == 0,
    optim_ha_nit = ha_results$counts,
    function_time = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  )
}

remove_highly_correlated <- function(ld_matrix, snp_ordering, max_correlation) {
  # Remove NAs from the correlation matrix
  ld_matrix[is.na(ld_matrix)] <- 0

  # Initialize a logical vector to mark rows/columns to remove
  idxs_to_remove <- rep(FALSE, nrow(ld_matrix))

  # Find pairs of SNPs that are highly correlated
  correlated_indices <- which(abs(ld_matrix) >= max_correlation, arr.ind = TRUE)

  to_remove <- vector("logical", length = nrow(ld_matrix))

  # Loop through correlated pairs and mark SNPs to remove
  for (i in seq_len(nrow(correlated_indices))) {
    a <- correlated_indices[i, 1]
    b <- correlated_indices[i, 2]

    if ((!to_remove[a] && !to_remove[b]) && (a != b)) {
      to_remove[b] <- TRUE
    }
  }

  # Remove highly correlated SNPs from both the matrix and the SNP list
  ld_matrix <- ld_matrix[!to_remove, !to_remove]
  pruned_snps <- snp_ordering[!to_remove]

  return(list(ld_matrix = ld_matrix, pruned_snps = pruned_snps))
}

mr_link2_analysis <- function(exposure_betas, outcome_betas, ld_matrix, n_exp, n_out, max_correlation = 0.99) {
  # Prune highly correlated SNPs
  snp_ordering <- seq_len(nrow(ld_matrix)) # Simple numbering for SNPs
  pruned_data <- remove_highly_correlated(ld_matrix, snp_ordering, max_correlation)
  pruned_ld_matrix <- pruned_data$ld_matrix
  pruned_snps <- pruned_data$pruned_snps

  # Eigenvalue decomposition
  eigen_decomp <- eigen(pruned_ld_matrix)
  eigenvalues <- eigen_decomp$values
  eigenvectors <- eigen_decomp$vectors

  # Calculate cumulative variance explained
  variance_explained <- cumsum(eigenvalues) / sum(eigenvalues)

  # Select the eigenvectors that explain at least 99% of the variance
  threshold_index <- min(which(variance_explained >= 0.99))
  selected_eigenvectors <- eigenvectors[, 1:threshold_index]
  selected_eigenvalues <- eigenvalues[1:threshold_index]

  # Call the existing MR-link-2 function
  # Assuming `mr_link2` is already defined in R and takes the following arguments:
  # - selected_eigenvalues
  # - selected_eigenvectors
  # - exposure_betas
  # - outcome_betas
  # - n_exp (number of individuals in exposure dataset)
  # - n_out (number of individuals in outcome dataset)
  sigma_exp_guess <- 0.01 # You can adjust or estimate this as needed
  sigma_out_guess <- 0.001 # You can adjust or estimate this as needed

  # Subset the betas based on pruned SNPs
  exposure_betas_pruned <- exposure_betas[pruned_snps]
  outcome_betas_pruned <- outcome_betas[pruned_snps]

  # Perform the MR-link-2 estimation
  mr_result <- mr_link2(
    selected_eigenvalues = selected_eigenvalues,
    selected_eigenvectors = selected_eigenvectors,
    exposure_betas = exposure_betas_pruned,
    outcome_betas = outcome_betas_pruned,
    n_exp = n_exp,
    n_out = n_out,
    sigma_exp_guess = sigma_exp_guess,
    sigma_out_guess = sigma_out_guess
  )

  return(mr_result)
}


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

ld_dataset <- arrow::open_dataset(args$ld_folder)
indik <- 1
for (eqtl_gene_indik in eqtl_genes) {
  if (isTRUE(args$mrlink2)) {
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
      R2 = NA,
      eQTL_lead_eQTL_Z = NA,
      eQTL_lead_GWAS_Z = NA,
      alpha = NA,
      `se(alpha)` = NA,
      `p(alpha)` = NA,
      sigma_y = NA,
      `se(sigma_y)` = NA,
      `p(sigma_y)` = NA,
      sigma_x = NA,
      alpha_h0_sigma_x = NA,
      alpha_h0_sigma_y = NA,
      alpha_h0_loglik = NA,
      sigma_y_h0_alpha = NA,
      sigma_y_h0_sigma_x = NA,
      sigma_y_h0_loglik = NA,
      ha_loglik = NA,
      optim_alpha_h0_success = NA,
      optim_alpha_h0_nit = NA,
      optim_sigma_y_h0_success = NA,
      optim_sigma_y_h0_nit = NA,
      optim_ha_success = NA,
      optim_ha_nit = NA,
      function_time = NA
    )[-1]
  } else {
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
      R2 = NA,
      eQTL_lead_eQTL_Z = NA,
      eQTL_lead_GWAS_Z = NA
    )[-1]
  }

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
    ref <- open_dataset(args$reference) %>%
      filter(variant_index %in% eqtl$variant_index) %>%
      select(variant_index, variant) %>%
      collect() %>%
      as.data.table()

    ref <- ref[order(variant_index), ]

    rownames(eqtl_mat) <- ref$variant

    gwas_loci_iterator <- 0
    length_gwas_loci <- nrow(inp_coloc_f)

    # Read in LD
    ld_dataset2 <- ld_dataset %>% filter(chr == as.integer(unique(eqtl$chromosome)))
    message("LD loaded!")
    ld_mat <- get_ld_matrix_wide(
      permuted_dataset = ld_dataset2,
      variants = as.integer(ref$variant_index)
    )
    message("LD read!")
    ld_mat <- ld_mat[order(rownames(ld_mat)), order(colnames(ld_mat))]

    ld_variant_indices <- colnames(ld_mat)

    colnames(ld_mat) <- ref$variant
    rownames(ld_mat) <- ref$variant

    if (!all(colnames(ld_mat) == rownames(eqtl_mat))) {
      message("Variants not aligned!")
      exit()
    }

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
        if (isTRUE(args$mrlink2)) {
          ld_mat_analysis <- ld_mat[rownames(ld_mat) %in% eqtl2$variant, colnames(ld_mat) %in% eqtl2$variant]

          if (!all(colnames(ld_mat_analysis) == eqtl2$variant)) {
            message("Variants not aligned!")
            exit()
          }
          # if(!all(colnames(ld_mat_analysis) == gwas2$variant_index)){message("Variants not aligned!"); exit()}
          message("Running MrLink2...")
          res_mrlink2 <- mr_link2_analysis(
            exposure_betas = eqtl2$beta,
            outcome_betas = gwas2$beta,
            ld_matrix = ld_mat_analysis,
            n_exp = median(eqtl2$sample_size),
            n_out = unique(gwas2$N),
            max_correlation = 0.99
          )
          message("Running MrLink2...done!")
          rm(ld_mat_analysis)
        }
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

                res_temp$R <- ld_mat[rownames(ld_mat) == res_temp$hit1, colnames(ld_mat) == res_temp$hit2]
                res_temp$eQTL_Z <- eqtl2[variant_index == ref[variant == res_temp$hit1]$variant_index]$beta / eqtl2[variant_index == ref[variant == res_temp$hit1]$variant_index]$standard_error
                res_temp$GWAS_Z <- gwas[variant_index == ref[variant == res_temp$hit1]$variant_index]$beta / gwas[variant_index == ref[variant == res_temp$hit1]$variant_index]$se
              } else {
                res_temp$coloc_cs_size <- NA
                res_temp$cs_variants <- NA
                res_temp$cs_PP <- NA
                res_temp$R <- NA
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
                res_temp$R <- ld_mat[rownames(ld_mat) == res_temp$hit1, colnames(ld_mat) == res_temp$hit2]
                res_temp$eQTL_Z <- eqtl2[variant_index == ref[variant == res_temp$hit1]$variant_index]$beta / eqtl2[variant_index == ref[variant == res_temp$hit1]$variant_index]$standard_error
                res_temp$GWAS_Z <- gwas[variant_index == ref[variant == res_temp$hit1]$variant_index]$beta / gwas[variant_index == ref[variant == res_temp$hit1]$variant_index]$se
              } else {
                res_temp$coloc_cs_size <- NA
                res_temp$cs_variants <- NA
                res_temp$cs_PP <- NA
                res_temp$R <- NA
                res_temp$eQTL_Z <- NA
                res_temp$GWAS_Z <- NA
              }
            }
            if (ncol(res_temp) == 21) {
              if (isTRUE(args$mrlink2)) {
                res_temp <- data.table(res_temp, as.data.table(res_mrlink2)[1, ])
              }

              if (isTRUE(args$mrlink2)) {
                res_temp <- res_temp[, c(1:8, 14:15, 9:13, 16:21, 22:ncol(res)), with = FALSE]
                colnames(res_temp) <- c(
                  "eQTL_gene", "eQTL_locus", "eQTL_type", "GWAS", "GWAS_locus", "nsnps", "hit1", "hit2",
                  "signal1", "signal2", "PP_H0", "PP_H1", "PP_H2", "PP_H3", "PP_H4", "coloc_cs95_size",
                  "cs95_variants", "cs95_PIP", "R2", "eQTL_lead_eQTL_Z", "eQTL_lead_GWAS_Z", colnames(res_temp)[22:ncol(res_temp)]
                )
              } else {
                res_temp <- res_temp[, c(1:8, 14:15, 9:13, 16:21), with = FALSE]
                colnames(res_temp) <- c(
                  "eQTL_gene", "eQTL_locus", "eQTL_type", "GWAS", "GWAS_locus", "nsnps", "hit1", "hit2",
                  "signal1", "signal2", "PP_H0", "PP_H1", "PP_H2", "PP_H3", "PP_H4", "coloc_cs95_size",
                  "cs95_variants", "cs95_PIP", "R2", "eQTL_lead_eQTL_Z", "eQTL_lead_GWAS_Z"
                )
              }
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
    rm(eqtl, eqtl_mat, ld_mat, ld_dataset2)
    gc()
    message(paste("Finalised with eQTL iterations"))
  }
}
warnings()
