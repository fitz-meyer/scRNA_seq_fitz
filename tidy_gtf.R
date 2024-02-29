library(tidyverse)

setwd("/Users/emilyfitzmeyer/Desktop/scRNAseq/analyses/cxt_annotation/all_CDS/")
col_names <- c("seqregion", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
gtf <- read.delim("cxt_genome_CDS_readInR.gtf", col.names = col_names)

sel_gtf <- gtf %>%
  select(1, 3, 4, 5, 9)

sep_gtf <- tidyr::separate(sel_gtf, col = attribute, into = c("A", "B", "C", "D", "E", "F", "G"), sep = ";" )

sel2_gtf <- sep_gtf %>%
  select(1, 2, 3, 4, 5)

sel2_gtf <- dplyr::rename(sel2_gtf, gene_id = A)

sel2_gtf <- sel2_gtf %>%
  mutate(gene_id = str_replace(sel2_gtf$gene_id, "gene_id ", ""))

write.csv(sel2_gtf, "tidy_cxt_CDS.csv", row.names = FALSE)
