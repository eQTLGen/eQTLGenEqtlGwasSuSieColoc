#!/usr/bin/env Rscript

library(data.table)
library(rtracklayer)
library(argparse)

setDTthreads(1)
parser <- ArgumentParser(description = 'Run coloc analysis.')


parser$add_argument('--gtf', 
  metavar = 'file', 
  type = 'character',
  help = 'GTF file.')

args <- parser$parse_args()

gtf <- readGFF(args$gtf)
gtf <- as.data.table(gtf)
gtf <- unique(gtf[type %in% "gene", c(9, 11), with = FALSE])

fwrite(gtf, "GeneNames.txt", sep = "\t")
