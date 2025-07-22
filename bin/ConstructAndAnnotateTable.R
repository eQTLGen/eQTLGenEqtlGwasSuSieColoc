#!/usr/bin/env Rscript

library(data.table)
library(rtracklayer)
library(stringr)
library(dplyr)
library(argparse)

parser <- ArgumentParser(description = "Make annotation table.")

parser$add_argument("--finemapped_eqtl",
  metavar = "file",
  type = "character",
  help = "Directory containing finemapped eQTL loci."
)

parser$add_argument("--finemapped_gwas",
  metavar = "file",
  type = "character",
  help = "Directory containing finemapped GWAS loci."
)

parser$add_argument("--genefilter",
  metavar = "file",
  type = "character",
  help = "File with ENSG IDs to confine the analysis to."
)

args <- parser$parse_args()

genefilter <- fread(args$genefilter, header = FALSE)

print(head(genefilter))
print(nrow(genefilter))

eqtl <- list.files(args$finemapped_eqtl)
gwas <- list.files(args$finemapped_gwas)

eqtl_file <- eqtl
eqtl_gene <- str_remove(eqtl, "__.*")
eqtl_type <- str_remove(eqtl, ".*___") %>% str_remove(".txt.gz")
eqtl_chr <- str_remove(eqtl, "___.*") %>% str_remove(".*__") %>% str_remove(":.*")
eqtl_start <- str_remove(eqtl, "___.*") %>% str_remove(".*__") %>% str_remove(".*:") %>% str_remove("_.*")
eqtl_end <- str_remove(eqtl, "___.*") %>% str_remove(".*__") %>% str_remove(".*:") %>% str_remove(".*_")

print(head(eqtl))
print(head(eqtl_gene))
print(head(eqtl_type))
print(head(eqtl_chr))
print(head(eqtl_start))
print(head(eqtl_end))

eqtl_table <- data.table(
    gene = eqtl_gene,
    type = eqtl_type,
    loc_chr = as.integer(eqtl_chr),
    loc_start = as.numeric(eqtl_start),
    loc_end = as.numeric(eqtl_end),
    file = eqtl_file
    )

print(head(eqtl_table, 3))

if (nrow(genefilter) > 0){eqtl_table <- eqtl_table[gene %in% genefilter$V1]}
message(paste(length(unique(eqtl_table$gene)), "genes in the analysis!"))

gwas_file <- gwas
gwas_gene <- str_remove(gwas, "__.*")
gwas_type <- str_remove(gwas, ".*___") %>% str_remove(".txt.gz")
gwas_chr <- str_remove(gwas, "___.*") %>% str_remove(".*__") %>% str_remove(":.*")
gwas_start <- str_remove(gwas, "___.*") %>% str_remove(".*__") %>% str_remove(".*:") %>% str_remove("_.*")
gwas_end <- str_remove(gwas, "___.*") %>% str_remove(".*__") %>% str_remove(".*:") %>% str_remove(".*_")

gwas_table <- data.table(
    gene = gwas_gene,
    type = gwas_type,
    loc_chr = as.integer(gwas_chr),
    loc_start = as.numeric(gwas_start),
    loc_end = as.numeric(gwas_end),
    file = gwas_file
    )

print(head(gwas_table, 40))

eqtl_table <- rbind(eqtl_table, gwas_table)

print(summary(eqtl_table))

# Annotate
eqtl_table_cis <- eqtl_table[type %in% c("cis", "trans")]
eqtl_table_trans <- eqtl_table[type == "gwas"]

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

colnames(overlaps) <- c("eqtl_gene", "loc_chr", "eqtl_loc_start", "eqtl_loc_end", "eqtl_file",
                        "gwas", "gwas_loc_start", "gwas_loc_end", "gwas_file")

message(paste(nrow(overlaps), "rows in the output table."))

fwrite(overlaps, "EqtlGwasOverlaps.txt", sep = "\t")
