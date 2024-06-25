library(scCustomize)
library(tidyverse)
library(Seurat)


merged.srt <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq_pub/mergeDPI_oldNorm/alldpi_oldNorm_nmn_res0.6.rds")
#rename clusters
new_idents <- c("EC-like-1", "1", "EC-like-2", "ISC/EB", "4", "EC", "EC-like-3", "VM-1", "8", "VM-2", "ISC/EB-prol", "cardia-1", "EE", "cardia-2", "HC-1", "MT", "HC-2", "17", "cardia-prol", "19")
merged.srt <- Rename_Clusters(merged.srt, new_idents = new_idents)

#add cell type info to metadata
ident_num <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19")
cell_type <- c("EC-like-1", "1", "EC-like-2", "ISC/EB", "4", "EC", "EC-like-3", "VM-1", "8", "VM-2", "ISC/EB-prol", "cardia-1", "EE", "cardia-2", "HC-1", "MT", "HC-2", "17", "cardia-prol", "19")
tib <- tibble(ident_num, cell_type)
merged.srt <- Add_Sample_Meta(merged.srt, meta_data = tib, join_by_seurat = "seurat_clusters", join_by_meta = "ident_num")

my_levels <- c("VM-2", "VM-1", "MT", "ISC/EB-prol", "ISC/EB", "HC-2", "HC-1", "EE", "EC-like-3", "EC-like-2", "EC-like-1", "EC", "cardia-2", "cardia-1", "cardia-prol", "19", "17", "8", "4", "1")
factor(Idents(merged.srt), levels = my_levels)
Idents(merged.srt) <- factor(Idents(merged.srt), levels = my_levels)

meta.data <- merged.srt[[]]

counts <- group_by(meta.data, condition, cell_type) %>% summarise(count = n())

ggplot(counts, aes(cell_type, count, fill = condition)) +
  theme(axis.text.x = element_text(size = 13, angle = 45, vjust = 0.5)) +
  theme(axis.text.y = element_text(size = 15)) +
  geom_bar(stat = 'identity')
