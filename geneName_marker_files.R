library(tidyverse)

setwd("/Users/emilyfitzmeyer/Desktop/temp/")
files <- Sys.glob("*.csv")
filenames <- str_replace(files, ".csv", "_named.csv")

list <- lapply(files, read_csv)
names(list) <- filenames

gene_name_list <- read_csv("/Users/emilyfitzmeyer/Desktop/scRNAseq/cxt_annotation/all_CDS_eggNOG/aggregated_geneID_name_list.csv")

merge_fun <- function(x) {
  merge(x, gene_name_list, by = "gene_ID", all.x = TRUE)
}

merged_list <- map(list, merge_fun)

setwd("/Users/emilyfitzmeyer/Desktop/")

for(i in 1:length(merged_list)){
  write.csv(merged_list[[i]], file = filenames[i], row.names = FALSE)
}




