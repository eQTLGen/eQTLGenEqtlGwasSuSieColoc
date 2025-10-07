#!/usr/bin/env Rscript

library(coloc)
library(susieR)
library(arrow)
library(data.table)
library(stringr)
library(dplyr)
library(argparse)
library(Rfast)
library(IGUtilityPackage)
library(duckdbfs)

set.seed(8)
options(scipen = 999)

setDTthreads(1)
parser <- ArgumentParser(description = "Run finamap analysis for GWAS.")

parser$add_argument("--gwas_info_file",
  metavar = "file",
  type = "character",
  help = "File containing GWAS parquet name, GWAS type (quant/cc), and effective N."
)

parser$add_argument("--maf_file",
  metavar = "file",
  type = "character",
  help = "File containing variant index and MAF."
)

parser$add_argument("--ld_folder",
  metavar = "file",
  type = "character",
  help = "Folder containing permuted LD files."
)

parser$add_argument("--full_lbf_matrix",
                    default = F,
                    action = 'store_true',
                    help="Return the entire LBF matrix instead of just those forming credible sets by providing this flag"
)

args <- parser$parse_args()

runsusie <- function(d,suffix=1,
                  maxit=100,repeat_until_convergence=TRUE,
                  s_init=NULL, ...) {
  coloc:::check_dataset(d,suffix,req=c("beta","varbeta","LD","snp","N"))
  coloc:::check_ld(d,d$LD)

  LD=d$LD[d$snp,d$snp] # just in case
  snp=d$snp

  converged=FALSE;
  ## set some defaults for susie arguments
  susie_args=list(...)
  if("max_iter" %in% names(susie_args)) {
    maxit=susie_args$max_iter
    susie_args = susie_args[ setdiff(names(susie_args), "max_iter") ]
  }
  ## at 0.12.6 susieR introduced need for n = sample size
  if(!("n" %in% names(susie_args)))
    susie_args=c(list(n=d$N), susie_args)

  while(!converged) {
    message("running max iterations: ",maxit)
    str(c(list(bhat = d$beta, shat = sqrt(d$varbeta), R=LD, max_iter=maxit), susie_args))
    res=do.call(susie_rss,
                c(list(bhat = d$beta, shat = sqrt(d$varbeta), R=LD, max_iter=maxit), susie_args))
    converged=res$converged; #s_init=res; maxit=maxit*2
    message("\tconverged: ",converged)
    if(!converged && repeat_until_convergence==FALSE)
      stop("susie_rss() did not converge in ",maxit," iterations. Try running with run_until_convergence=TRUE")
    if(!converged)
      maxit=maxit * 100 # no point in half measures!
  }
  res=annotate_susie(res, snp, LD)

  return(res)
}

# functions
safe_runsusie <- function(...) {
  tryCatch(
    {
      result <- runsusie(...)
      if (is.null(result) || isFALSE(result$converged)) {
        message("runsusie did not converge.")
        return(NULL)
      }
      return(result)
    },
    error = function(e) {
      message("runsusie failed with error: ", e$message)
      return(NULL)
    }
  )
}

compute_yty <- function(beta, se, p, R, n, k) {
  # beta and se should be standarised
  beta_s <- beta * sqrt(2 * p * (1 - p))
  se_s <- se * sqrt(2 * p * (1 - p))

  # Y'Y =  Bj^2 (Xj'Xj) + Var(Bj)(Xj'Xj)(N - k)
  XjtXj <- (n - 1) * diag(R)
  yty <- beta_s**2 * XjtXj + se_s**2 * XjtXj * (n - k)

  return(median(yty))
}

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

print(args)

locus_files <- list.files(pattern = "___temp.txt.gz")
print(locus_files)
gwas_info <- fread(args$gwas_info_file)

message("Reading MAF file...")
maf <- fread(args$maf)
maf <- maf[order(variant_index)]
message("Reading MAF file...done!")

# LD matrix
message("Opening LD file...")
ld_dataset <- duckdbfs::open_dataset(args$ld_folder)
message("Opening LD file...done!")

filter_lbf_mat <- !args$full_lbf_matrix
if (filter_lbf_mat) {
  message("Filtering LBF matrix for those forming credible sets!")
} else {
  message("Keeping the full LBF matrix without filtering for those forming credible sets!")
}

for (i in 1:length(locus_files)) {
  if (exists("gwas_susie")) {
    rm(gwas_susie)
  }
  if (exists("lbf")) {
    rm(lbf)
  }
  if (exists("trace_gwas")) {
    rm(trace_gwas)
  }

  remove_maf <- 0

  locus <- fread(locus_files[i])
  pheno_id <- str_replace(locus_files[i], "__.*", "")

  locus_id <- str_replace(locus_files[i], "___temp.*", "")
  print(locus_id)
  locus_id <- str_replace(locus_id, ".*__", "")

  gwas_type <- gwas_info[pheno == pheno_id]$type

  message(paste0("Analysing: ", locus_id, " ", i, "/", length(locus_files), "..."))

  start.time <- Sys.time()

  message(paste0("Removing variants with Beta == 0.0 ", nrow(locus[beta == 0])))
  locus <- locus[beta != 0]

  temp_chr <- as.integer(unique(locus$chromosome))
  # temp_start <- min(as.integer(locus$bp))
  # temp_end <- max(as.integer(locus$bp))

  message(paste(nrow(locus), "variants"))

  ld_dataset2 <- ld_dataset %>% filter(chr == temp_chr)

  ld_mat <- get_ld_matrix_wide(
    permuted_dataset = ld_dataset2,
    variants = as.integer(locus$variant_index)
  )

  print(dim(ld_mat))

  if (ncol(ld_mat) > 0) {

    locus$variant_index <- as.character(locus$variant_index)
    locus <- locus[variant_index %in% rownames(ld_mat)]

    ld_mat <- ld_mat[order(as.integer(rownames(ld_mat))), order(as.integer(colnames(ld_mat))), drop=F]
    locus <- locus[order(as.integer(variant_index))]

    # Remove sumstat variants that are not in LD matrix.
    locus <- locus[as.integer(locus$variant_index) %in% as.integer(colnames(ld_mat))]
    locus <- locus[!duplicated(locus$variant_index)]

    print(locus[duplicated(variant_index)])
    print(dim(ld_mat))

    if (ncol(ld_mat) == 0) {
      message("Fine-mapping did not succeed. Locus is empty (after filtering)!")
      next
    }

    # Compare Zs 
    cond_z <- kriging_rss(locus$beta / locus$se, ld_mat, n = median(locus$N))
    cond_z <- cond_z$conditional_dist

    Ps <- ZtoP(cond_z$z_std_diff)
    message(paste0("Removing variants with abs Z-score difference of >", min(abs(cond_z[Ps <= (0.05 / nrow(locus)), ]$z_std_diff))))
    message(paste0(nrow(locus[Ps <= (0.05 / nrow(locus))]), "/", nrow(locus), " variants removed due to outlier Z-score."))

    remove_outlier <- nrow(locus[Ps <= (0.05 / nrow(locus))])


    ld_mat <- ld_mat[Ps > (0.05 / nrow(locus)), Ps > (0.05 / nrow(locus)), drop=F]
    locus <- locus[Ps > (0.05 / nrow(locus))]

    if (ncol(ld_mat) == 0) {
      message("Fine-mapping did not succeed. Locus is empty (after filtering)!")
      next
    }

    maf2 <- maf[as.character(variant_index) %in% as.character(colnames(ld_mat))]

    print(nrow(locus))
    print(dim(ld_mat))
    print(nrow(maf2))

    print(locus[!variant_index %in% colnames(ld_mat)])
    print(maf2[!variant_index %in% colnames(ld_mat)])

    ld_mat <- ld_mat[order(as.integer(rownames(ld_mat))), order(as.integer(colnames(ld_mat))), drop=F]
    locus <- locus[order(as.integer(variant_index))]
    maf2 <- maf2[order(as.integer(variant_index))]

    print(locus[variant_index != colnames(ld_mat)])
    print(maf2[variant_index != colnames(ld_mat)])
    print(locus[variant_index != maf2$variant_index])

    locus$MAF <- maf2$MAF

    if(!all(colnames(ld_mat) == locus$variant_index)){stop("ERROR: Column names of LD matrix do not match locus variant indices!")}
    if(!all(as.character(colnames(ld_mat)) == as.character(maf2$variant_index))){stop("ERROR: Column names of LD matrix do not match MAF variant indices!")}

    Neff <- gwas_info[pheno == pheno_id]$Neff

    # SuSiE analysis:
    if (gwas_type == "quant"){

    standardised_pheno <- gwas_info[pheno == pheno_id]$standardised

    if (standardised_pheno == "yes"){

    coloc_inp_gwas <- list(
      beta = locus$beta,
      varbeta = locus$se * locus$se,
      snp = as.character(locus$variant_index),
      position = as.integer(locus$bp),
      N = Neff,
      type = gwas_type,
      sdY = 1,
      LD = ld_mat
    )
    } else {

      message(paste("Removing due to MAF=0", nrow(locus[MAF == 0 | MAF == 1]), "variants"))
      locus <- locus[MAF > 0 & MAF < 1]

      remove_maf <- nrow(locus[MAF == 0 | MAF == 1])

      ld_mat <- ld_mat[as.character(rownames(ld_mat)) %in% as.character(locus$variant_index), as.character(colnames(ld_mat)) %in% as.character(locus$variant_index), drop=F]

      if(!all(colnames(ld_mat) == locus$variant_index)){stop("ERROR: Column names of LD matrix do not match locus variant indices!")}

      coloc_inp_gwas <- list(
      beta = locus$beta,
      varbeta = locus$se * locus$se,
      snp = as.character(locus$variant_index),
      position = as.integer(locus$bp),
      N = Neff,
      MAF = as.numeric(locus$MAF),
      type = gwas_type,
      LD = ld_mat
      )
    }
    } else if(gwas_type == "cc"){

      Neff <- gwas_info[pheno == pheno_id]$Neff

      coloc_inp_gwas <- list(
        beta = locus$beta,
        varbeta = locus$se * locus$se,
        snp = as.character(locus$variant_index),
        position = as.integer(locus$bp),
        N = Neff,
        MAF = as.numeric(locus$MAF),
        type = gwas_type,
        LD = ld_mat
      )
    }

    # Calculate lambda
    lambda <- estimate_s_rss(locus$beta / locus$se, ld_mat, n = median(locus$N))
    message(paste("Lambda:", lambda))

    var_y <- coloc:::sdY.est(locus$se^2, locus$MAF, median(locus$N))

    gwas_susie <- safe_runsusie(coloc_inp_gwas, var_y = var_y, L = 10, estimate_residual_variance = FALSE, repeat_until_convergence = FALSE, maxit = 10000)
    if (!is.null(gwas_susie)) {
      trace_gwas <- list(L = 10, estimate_residual_variance = FALSE, niter = gwas_susie$niter)
    }
    if (is.null(gwas_susie) | isFALSE(gwas_susie$converged) | is.null(gwas_susie$sets$cs)) {
      message("None of the nCS values led to convergence. Retrying with different L.")
      message("Attempting to run SuSie with L = ", 5, ", and estimate_residual_variance = FALSE")
      gwas_susie <- safe_runsusie(coloc_inp_gwas, var_y = var_y, L = 5, estimate_residual_variance = FALSE, repeat_until_convergence = FALSE, maxit = 10000)
      if (!is.null(gwas_susie)) {
        trace_gwas <- list(L = 5, estimate_residual_variance = FALSE, niter = gwas_susie$niter)
      }
    }
    if (is.null(gwas_susie) | isFALSE(gwas_susie$converged) | is.null(gwas_susie$sets$cs)) {
      message("None of the nCS values led to convergence. Retrying with different L.")
      message("Attempting to run SuSie with L = ", 3, ", and estimate_residual_variance = FALSE")
      gwas_susie <- safe_runsusie(coloc_inp_gwas, var_y = var_y, L = 3, estimate_residual_variance = FALSE, repeat_until_convergence = FALSE, maxit = 10000)
      if (!is.null(gwas_susie)) {
        trace_gwas <- list(L = 3, estimate_residual_variance = FALSE, niter = gwas_susie$niter)
      }
    }

    if (!is.null(gwas_susie) & !isFALSE(gwas_susie$converged) & !is.null(gwas_susie$sets$cs)) {
      # Get LBF matrix including all credible sets
      lbf_matrix <- t(gwas_susie$lbf_variable)
      colnames(lbf_matrix) <- paste0("lbf_cs_", 1:(ncol(lbf_matrix)))

      str(lbf_matrix)

      # Identify the credible sets that are identified by SuSiE (after internal deduplication steps)
      credible_set_names <- names(gwas_susie$sets$cs)
      credible_set_ids <- str_replace(credible_set_names, "L", "lbf_cs_")

      lbf_matrix_filtered <- lbf_matrix[, c(credible_set_ids), drop = F]

      # Construct final dataframe with results, including LBF matrix
      if (filter_lbf_mat) {
        lbf_matrix <- lbf_matrix_filtered
      }

      str(lbf_matrix)
      print(colnames(lbf_matrix))

      lbf <- data.table(
        variant_index = colnames(gwas_susie$lbf_variable),
        L = trace_gwas$L,
        niter = trace_gwas$niter,
        beta = locus$beta,
        se = locus$se,
        MAF = locus$MAF,
        N = Neff, lambda = lambda,
        lbf_matrix[, colnames(lbf_matrix), drop=F])

      # Add credible set identifiers to the output file, as well as a column with filtered lbf identifiers
      lbf$pip <- gwas_susie$pip
      lbf$CS <- NA
      lbf$CS_pass <- NA
      if(length(gwas_susie$sets$cs) != 0){
        for(l in 1:length(gwas_susie$sets$cs)){
          credible_set_id_as_numeric <- gsub("L", "", credible_set_names[l])
          lbf$CS[gwas_susie$sets$cs[[l]]] <- credible_set_id_as_numeric
          if (max(lbf_matrix_filtered[,credible_set_ids[l]]) > 2) {
            lbf$CS_pass[gwas_susie$sets$cs[[l]]] <- credible_set_id_as_numeric
          }
        }
      }

      fwrite(lbf, paste0(pheno_id, "__", locus_id, "___", "gwas", ".txt.gz"), sep = "\t", row.names=F)

      fwrite(
        data.table(
          gene = pheno_id,
          type = "gwas",
          loc_chr = unique(locus$chromosome),
          loc_start = min(locus$bp),
          loc_end = max(locus$bp),
          remove_outlier = remove_outlier,
          remove_maf = remove_maf,
          file = paste0(pheno_id, "__", locus_id, "___", "gwas", ".txt.gz")
        ),
        paste0(pheno_id, "__", locus_id, "___", "gwas", "_AnnotationFile.txt"),
        sep = "\t"
      )
    } else {

      gwas_finemap_abf <- as.data.frame(finemap.abf(coloc_inp_gwas))
      rownames(gwas_finemap_abf) <- gwas_finemap_abf$snp

      # Construct the SuSiE-style LBF table
      lbf <- data.table(
        variant_index = locus$variant_index,
        L = NA,
        niter = NA,
        beta = locus$beta,
        se = locus$se,
        MAF = locus$MAF,
        N = Neff,
        lambda = lambda,
        lbf_cs_1 = gwas_finemap_abf[as.character(locus$variant_index), ]$lABF.
      )

      # Add credible set identifiers to the output file, as well as a column with filtered lbf identifiers
      lbf$pip <- gwas_finemap_abf[as.character(locus$variant_index), ]$SNP.PP
      lbf$CS <- NA
      lbf$CS_pass <- NA

      # Sort by PIP descending and compute cumulative sum
      setorder(lbf, -pip)
      lbf[, cumulative_pip := cumsum(pip)]

      # Select variants until cumulative_PIP reaches threshold
      threshold <- 0.95
      variants_passed_threshold <- which(lbf[, cumulative_pip > threshold])

      if (length(variants_passed_threshold) == 0) {
        message("Total probability not greater than threshold")
      } else {
        first_variant_passed_threshold <- variants_passed_threshold[1]

        lbf$CS[1:first_variant_passed_threshold] <- "1"
        if (max(lbf$lbf_cs_1) > 2) {
          lbf$CS_pass[1:first_variant_passed_threshold] <- "1"
        }

        lbf <- lbf[order(as.numeric(lbf$variant_index)), ]
        lbf$cumulative_pip <- NULL

      }


      fwrite(lbf, paste0(pheno_id, "__", locus_id, "___", "gwas", ".txt.gz"), sep = "\t")

      fwrite(
        data.table(
          gene = pheno_id,
          type = "gwas",
          loc_chr = unique(locus$chromosome),
          loc_start = min(locus$bp),
          loc_end = max(locus$bp),
          remove_outlier = remove_outlier,
          remove_maf = remove_maf,
          file = paste0(pheno_id, "__", locus_id, "___", "gwas", ".txt.gz")
        ),
        paste0(pheno_id, "__", locus_id, "___", "gwas", "_AnnotationFile.txt"),
        sep = "\t"
      )

      message("SuSiE fine-mapping did not succeed, performed Finemap ABF")
    }
  } else {
   
    message("Fine-mapping did not succeed due to missing LD matrix.")
  }
  end.time <- Sys.time()
  message(paste("Elapsed:", end.time - start.time))
  message(paste0("Analysing: ", locus_id, " ", i, "/", length(locus_files), "...done!"))
  rm(ld_mat, maf2, locus, gwas_susie, ld_dataset2, cond_z)
  gc()
}

warnings()
