library(tidyverse)

setwd("/Users/emilyfitzmeyer/Desktop/scRNAseq/cxt_annotation/all_CDS_eggNOG/")
gene_id_df <- read.csv("tidy_cxt_CDS.csv")
eggnog_df <- read.csv("eggnog_feature_tidy_cxt_CDS.csv")

gene_id_df <- gene_id_df %>% 
  mutate(source = paste(gene_id_df$seqregion, gene_id_df$feature, gene_id_df$start, gene_id_df$end, sep = "_")) %>%
  select(6, 5)

eggnog_df <- eggnog_df %>%
  mutate(source = paste(eggnog_df$seqregion, eggnog_df$feature, eggnog_df$start, eggnog_df$end, sep = "_")) %>%
  select(25, 5, 10, 11, 12, 24)

df_merge <- merge(gene_id_df, eggnog_df, by = "source", all.x = TRUE, all.y = FALSE)

write.csv(df_merge, "geneID_eggnog_merge_WITHSEED.csv", row.names = FALSE)

