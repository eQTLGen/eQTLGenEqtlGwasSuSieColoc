#!/usr/bin/env Rscript

library(data.table)
library(rtracklayer)
library(stringr)
library(dplyr)
library(argparse)

parser <- ArgumentParser(description = "Make annotation table.")

parser$add_argument(
  "--finemapped_gwas",
  metavar = "file",
  type = "character",
  help = "Directory containing finemapped eQTL loci."
)

parser$add_argument(
  "--a", nargs="+",
  metavar = "character",
  type = "character",
  help = "GWAS traits to test"
)

parser$add_argument(
  "--b", nargs="+",
  metavar = "character",
  type = "character",
  help = "GWAS traits to test."
)


args <- parser$parse_args()

gwas_filter_a <- args$a
gwas_filter_b <- args$b

print(head(gwas_filter_a))
print(nrow(gwas_filter_b))

gwas <- list.files(args$finemapped_gwas)

gwas_file <- gwas
gwas_gene <- str_remove(gwas, "__.*")
gwas_type <- str_remove(gwas, ".*___") %>% str_remove(".txt.gz")
gwas_chr <- str_remove(gwas, "___.*") %>%
  str_remove(".*__") %>%
  str_remove(":.*")
gwas_start <- str_remove(gwas, "___.*") %>%
  str_remove(".*__") %>%
  str_remove(".*:") %>%
  str_remove("_.*")
gwas_end <- str_remove(gwas, "___.*") %>%
  str_remove(".*__") %>%
  str_remove(".*:") %>%
  str_remove(".*_")

gwas_table <- data.table(
  gene = gwas_gene,
  type = gwas_type,
  loc_chr = as.integer(gwas_chr),
  loc_start = as.numeric(gwas_start),
  loc_end = as.numeric(gwas_end),
  file = gwas_file
)

print(gwas_table %>% filter(is.na(loc_chr) | is.na(loc_start) | is.na(loc_end)))

print(head(gwas_table, 40))

print(summary(gwas_table))

# Annotate
eqtl_table_cis <- gwas_table[gene %in% gwas_filter_a]
eqtl_table_trans <- gwas_table[gene %in% gwas_filter_b]

eqtl_table_cis[, loc_start := as.integer(loc_start)]
eqtl_table_cis[, loc_end := as.integer(loc_end)]
eqtl_table_cis[, loc_chr := as.integer(loc_chr)]

eqtl_table_trans[, loc_start := as.integer(loc_start)]
eqtl_table_trans[, loc_end := as.integer(loc_end)]
eqtl_table_trans[, loc_chr := as.integer(loc_chr)]

setkey(eqtl_table_cis, loc_chr, loc_start, loc_end)
setkey(eqtl_table_trans, loc_chr, loc_start, loc_end)

overlaps <- foverlaps(eqtl_table_trans, eqtl_table_cis, nomatch = 0L)

overlaps <- overlaps[, c(2, 1, 4:7, 9:11), with = FALSE]

colnames(overlaps) <- c("gwas_a", "loc_chr", "a_loc_start", "a_loc_end", "a_file",
                        "gwas_b", "b_loc_start", "b_loc_end", "b_file")

overlaps <- overlaps[order(gwas_a)]

message(paste(nrow(overlaps), "rows in the output table."))

fwrite(overlaps, "GwasPairsOverlaps.txt", sep = "\t")
