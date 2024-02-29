library(scCustomize)
library(tidyverse)
library(Seurat)

srt <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/merged_rds/4dpi_res0.6.rds")
split.srt <- SplitObject(srt, split.by = "sample")

AverageExpression(split.srt$wnvMg7, group.by = "ident", slot = 'data', features = "nbis-gene-2-utr")

Percent_Expressing(split.srt$wnvMg8, features = "nbis-gene-2-utr")

prop.table(table(Idents(split.srt$wnvMg9)))

cluster_proportions <- prop.table(table(Idents(split.srt$mg4c)))
as.data.frame(cluster_proportions)




wnvmg6 <- split.srt$wnvMg6
cl12_wnvmg6.srt <- subset(wnvmg6, idents = "12")
cl12_vRNA_mg6 <- FetchData(cl12_wnvmg6.srt, vars = "nbis-gene-2-utr")

cl13_wnvmg6.srt <- subset(wnvmg6, idents = "13")
cl13_vRNA_mg6 <- FetchData(cl13_wnvmg6.srt, vars = "nbis-gene-2-utr")

cl10_wnvmg6.srt <- subset(wnvmg6, idents = "10")
cl10_vRNA_mg6 <- FetchData(cl10_wnvmg6.srt, vars = "nbis-gene-2-utr")





wnvmg7 <- split.srt$wnvMg7
cl12_wnvmg7.srt <- subset(wnvmg7, idents = "12")
cl12_vRNA_mg7 <- FetchData(cl12_wnvmg7.srt, vars = "nbis-gene-2-utr")

cl13_wnvmg7.srt <- subset(wnvmg7, idents = "13")
cl13_vRNA_mg7 <- FetchData(cl13_wnvmg7.srt, vars = "nbis-gene-2-utr")

cl10_wnvmg7.srt <- subset(wnvmg7, idents = "10")
cl10_vRNA_mg7 <- FetchData(cl10_wnvmg7.srt, vars = "nbis-gene-2-utr")





wnvmg8 <- split.srt$wnvMg8
cl12_wnvmg8.srt <- subset(wnvmg8, idents = "12")
cl12_vRNA_mg8 <- FetchData(cl12_wnvmg8.srt, vars = "nbis-gene-2-utr")

cl13_wnvmg8.srt <- subset(wnvmg8, idents = "13")
cl13_vRNA_mg8 <- FetchData(cl13_wnvmg8.srt, vars = "nbis-gene-2-utr")

cl10_wnvmg8.srt <- subset(wnvmg8, idents = "10")
cl10_vRNA_mg8 <- FetchData(cl10_wnvmg8.srt, vars = "nbis-gene-2-utr")





wnvmg9 <- split.srt$wnvMg9
cl12_wnvmg9.srt <- subset(wnvmg9, idents = "12")
cl12_vRNA_mg9 <- FetchData(cl12_wnvmg9.srt, vars = "nbis-gene-2-utr")

cl13_wnvmg9.srt <- subset(wnvmg9, idents = "13")
cl13_vRNA_mg9 <- FetchData(cl13_wnvmg9.srt, vars = "nbis-gene-2-utr")

cl10_wnvmg9.srt <- subset(wnvmg9, idents = "10")
cl10_vRNA_mg9 <- FetchData(cl10_wnvmg9.srt, vars = "nbis-gene-2-utr")








cl13_4dpi.srt <- subset(srt, idents = "13")
cl12_4dpi.srt <- subset(srt, idents = "12")

cec6016 <- VlnPlot(cl13_4dpi.srt, features = "gene6016", group.by = "condition", pt.size = 0.5) +
  theme(legend.position = "none") +
  ggtitle("Cecropin_")

cec6015 <- VlnPlot(cl13_4dpi.srt, features = "gene6015", group.by = "condition", pt.size = 0.5) +
  theme(legend.position = "none") +
  ggtitle("Cecropin")

cec10089 <- VlnPlot(cl13_4dpi.srt, features = "gene10089", group.by = "condition", pt.size = 0.5) +
  theme(legend.position = "none") +
  ggtitle("Cecropin")

plot <- cec6016 + cec6015 + cec10089 + patchwork::plot_layout(ncol = 3)

ggsave("temp.png", plot = plot, device = png(), scale = 1, width = 9, height = 5, dpi = 300)
dev.off()


#srt <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/dimreduc/mg9_c_dimreduc.rds")

immn.gene.pctExp <- Percent_Expressing(split.srt$wnvMg4, assay = 'RNA',
                                       features = c("gene10409", "gene9527", "gene3946"), entire_object = TRUE)

immn.gene.avgExp <- AverageExpression(split.srt$wnvMg4, assays = 'RNA', group.by = "orig.ident",
                                       features = c("gene10409", "gene9527", "gene3946"))

immn.gene.pctExp <- Percent_Expressing(split.srt$mg9c, assay = 'RNA', 
                                       features = c("gene1501", "gene7624", "gene10322", "gene11340", "gene3044",
                                                    "gene6258", "gene8444", "gene11364", "gene6091", "gene6092",
                                                    "gene5869", "gene5870", "gene10292", "gene10293", "gene3397",
                                                    "gene14111", "gene2893", "gene3913", "gene2450", "gene5218", 
                                                    "gene13590", "gene3988", "gene3989", "gene4535", "gene14091"),
                                       entire_object = TRUE)


immn.gene.avgExp <- AverageExpression(srt, assays = 'RNA', group.by = "orig.ident",
                                      features = c("gene1501", "gene7624", "gene10322", "gene11340", "gene3044",
                                                   "gene6258", "gene8444", "gene11364", "gene6091", "gene6092",
                                                   "gene5869", "gene5870", "gene10292", "gene10293", "gene3397",
                                                   "gene14111", "gene2893", "gene3913", "gene2450", "gene5218", 
                                                   "gene13590", "gene3988", "gene3989", "gene4535", "gene14091"))

#ML domain
ml_dom_pctExp <- Percent_Expressing(srt, assay = 'RNA', features = "gene7411", entire_object = TRUE)
ml_dom_avgExp <- AverageExpression(srt, assays = 'RNA', features = "gene7411", group.by = "orig.ident")

#WNV 5' UTR expressing values:

srt <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/dimreduc/wnv_mg9_dimreduc.rds")

wnv.pctExp <- Percent_Expressing(srt, assay = 'RNA', features = "nbis-gene-2-utr", entire_object = TRUE)
wnv.avgExp <- AverageExpression(srt, assays = 'RNA', features = "nbis-gene-2-utr", group.by = "ident")

