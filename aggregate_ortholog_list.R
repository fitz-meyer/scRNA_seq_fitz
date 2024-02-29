library(tidyverse)

setwd("/Users/emilyfitzmeyer/Desktop/scRNAseq/cxt_annotation/all_CDS_eggNOG/")
name_df <- read.csv("geneID_eggnog_merge_WITHSEED.csv")

agg_name <- aggregate(name_df[c('seed_ortholog', 'COG_category', 'Description', 'Preferred_name', 'PFAMs')], 
                      by = name_df['gene_id'], paste, collapse = ", ")

write_csv(agg_name, file = "/Users/emilyfitzmeyer/Desktop/aggregated_geneID_name_list_WITHSEED.csv")

