#!/usr/bin/env Rscript

library(data.table)
library(arrow)
library(dplyr)
library(rtracklayer)
library(argparse)

setDTthreads(1)
parser <- ArgumentParser(description = 'Overlap cis and trans loci.')

parser$add_argument('--input_file', metavar = 'file', type = 'character',
                    help = 'Input file that needs to be parsed.')
parser$add_argument('--reference', metavar = 'file', type = 'character',
                    help = 'eQTLGen SNP reference file in parquet format.')
parser$add_argument('--gtf', metavar = 'file', type = 'character',
                    help = "ENSEMBL .gtf file, needs to be hg38.")

args <- parser$parse_args()

and <- fread(args$input_file)

ref <- open_dataset(args$reference) %>%
  filter(variant_index %in% and$variant_index) %>%
  select(variant_index, chromosome, bp) %>%
  collect() %>%
  as.data.table()

and <- merge(and, ref, by = "variant_index")

and <- and %>%
  group_by(phenotype, locus_chromosome, locus_start) %>%
  mutate(locus_start = min(bp), locus_end = max(bp)) %>%
  group_by(phenotype, locus_chromosome, locus_start, SusieRss_CS) %>%
  mutate(lead_signal = case_when(!is.na(SusieRss_CS) & SusieRss_pip == max(SusieRss_pip) ~ "yes", .default = "no")) %>%
  as.data.table()

gtf <- readGFF(args$gtf)
gtf <- as.data.table(gtf)
gtf <- unique(gtf[type %in% "gene" & gene_id %in% and$phenotype, c(9, 1, 4, 5), with = FALSE])

and <- merge(and, gtf, by.x = "phenotype", by.y = "gene_id")

and <- and %>%
  group_by(phenotype, locus_chromosome, locus_start) %>%
  mutate(type = case_when(lead_signal == "yes" & chromosome == seqid & (abs(start - bp) < 1000000 | abs(end - bp) < 1000000) ~ "cis",
                          lead_signal == "yes" & (chromosome != seqid | (abs(start - bp) > 5000000 & abs(end - bp) > 5000000)) ~ "trans",
                          lead_signal == "yes" & chromosome == seqid & ((abs(start - bp) < 5000000 & abs(start - bp) > 1000000 ) | abs(end - bp) < 5000000 & abs(end - bp) > 1000000) ~ "interim",
                          .default = "NA")) %>%
  mutate(locus_type = case_when(any(type == "cis") ~ "cis", !any(type == "cis") ~ "trans")) %>%
  as.data.table()

for (gene in unique(and$phenotype)){

  and_p <- and[phenotype %in% gene]

  and_p$locus_id <- paste0(and_p$locus_chromosome, ":", and_p$locus_start, "_", and_p$locus_end)

  for (locus in unique(and_p$locus_id)){
    and_p2 <- and_p[locus_id %in% locus]
    type <- unique(and_p2$locus_type)
    fwrite(and_p2[, -ncol(and_p2), with = FALSE], paste0(gene, "__", locus, "___", type, ".txt.gz"), sep = "\t")

    fwrite(
      data.table(gene = gene, 
      type = type, 
      loc_chr = unique(and_p2$locus_chromosome), 
      loc_start = unique(and_p2$locus_start), 
      loc_end = unique(and_p2$locus_end), 
      file = paste0(gene, "__", locus, "___", type, ".txt.gz")),
    paste0(gene, "__", locus, "___", type, "_AnnotationFile.txt"), sep = "\t"
    )

  }

}
