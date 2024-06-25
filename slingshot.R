#trajectory analysis w/ slingshot
#Seurat workflow

library(singleCellTK)
library(tidyverse)
library(slingshot)
library(Seurat)
library(igraph)
library(tradeSeq)
library(gridExtra)
library(viridis)

sessionInfo()
#define color palette
pal <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"), RColorBrewer::brewer.pal(7, "Set3"))

#Read in pre-analyzed (guided clustering workflow up through UMAP reduction) Seurat object
merged.srt <- readRDS("rds/alldpi_oldNorm_nmn_res0.6.rds")
#merged.sce <- as.SingleCellExperiment(merged.srt)
DimPlot(merged.srt, group.by = "RNA_snn_res.0.6")

#Save seurat object components as separate matrices for slingshot
dimred <- merged.srt@reductions$umap@cell.embeddings
clustering <- merged.srt$RNA_snn_res.0.6
counts <- as.matrix(merged.srt@assays$RNA@counts[merged.srt@assays$RNA@var.features, ])

# #run default slingshot lineage identification
# set.seed(1)
# lineages <- getLineages(dimred, clusterLabels = clustering, omega = TRUE)
# lineages.dat <- as.SlingshotDataSet(lineages)
# 
# #plot the lineages
# par(mfrow = c(1, 2))
# plot(dimred[, 1:2], col = pal[clustering], cex = 0.5, pch = 20)
# for (i in levels(clustering)) {
#   text(mean(dimred[clustering == i, 1]), mean(dimred[clustering == i, 2]), labels = i, font = 2)
# }
# plot(dimred[, 1:2], col = pal[clustering], cex = 0.5, pch = 20)
# lines(lineages.dat, lwd = 3, col = "black")
# 
# #define starting cluster (defining ISC/EBs cluster 3)
# set.seed(1)
# lineages <- getLineages(dimred, clusterLabels = clustering, start.clus = "3", omega = TRUE, omega_scale = 1.4)
# lineages.dat <- as.SlingshotDataSet(lineages)
# 
# #plot
# par(mfrow=c(1,2))
# plot(dimred[,1:2], col = pal[clustering],  cex=.5, pch = 16)
# for(i in levels(clustering)){ 
#   text( mean(dimred[clustering==i,1]),
#         mean(dimred[clustering==i,2]), labels = i,font = 2) }
# plot(dimred, col = pal[clustering],  cex = .5, pch = 16)
# lines(lineages.dat, lwd = 3, col = 'black')
# 
# curves <- getCurves(lineages, thresh = 0.01, stretch = 0.8, allow.breaks = TRUE, shrink = TRUE)
# curves.dat <- as.SlingshotDataSet(curves)
# 
# plot(dimred, col = pal[clustering], cex = .5, asp = 1, pch = 16)
# lines(curves.dat, lwd = 3, col = "black")

#Test with just epithelial cells (EC, EE, cardia, ISC/EB)
merged.srt <- readRDS("Desktop/scRNAseq_pub/mergeDPI_oldNorm/alldpi_oldNorm_nmn_res0.6.rds")
epithelial.srt <- subset(merged.srt, idents = c("3", "10", "6", "2", "0", "5", "12", "13", "11", "18"))
#epithelial.sce <- as.SingleCellExperiment(epithelial.srt)
DimPlot(epithelial.srt, group.by = "RNA_snn_res.0.6")

dimred <- epithelial.srt@reductions$umap@cell.embeddings
clustering <- epithelial.srt$RNA_snn_res.0.6
counts <- as.matrix(epithelial.srt@assays$RNA@counts[epithelial.srt@assays$RNA@var.features, ])

#pull GOIs (cell type markers)
all_counts <- epithelial.srt@assays$RNA@counts
goi_sel <- c("gene13104", "gene2167", "gene11957", "gene10632")
goi_counts <- all_counts[goi_sel,]
goi_counts <- as.matrix(goi_counts)
#append to counts matrix
counts <- rbind(counts, goi_counts)

#define starting cluster (defining ISC/EBs cluster 3)
set.seed(1)
lineages <- getLineages(dimred, clusterLabels = clustering, start.clus = c("3", "10"), end.clus = c("12", "5"))
lineages.dat <- as.SlingshotDataSet(lineages)

#plot
par(mfrow=c(1,2))
plot(dimred[,1:2], col = pal[clustering],  cex=.5, pch = 16)
for(i in levels(clustering)){ 
  text( mean(dimred[clustering==i,1]),
        mean(dimred[clustering==i,2]), labels = i,font = 2) }
plot(dimred, col = pal[clustering],  cex = .5, pch = 16)
lines(lineages.dat, lwd = 3, col = 'black')

curves <- getCurves(lineages, thresh = 0.01, stretch = 0.8, allow.breaks = TRUE, shrink = TRUE)
curves.dat <- as.SlingshotDataSet(curves)

plot(dimred, col = pal[clustering], cex = .5, asp = 1, pch = 16)
lines(curves.dat, lwd = 3, col = "black")

#remove genes that aren't likely to be impactful to speed up computation
filt_counts <- counts[rowSums(counts > 5) > ncol(counts)/100, ]
dim(filt_counts)
#append to filt_counts matrix
filt_counts <- rbind(filt_counts, goi_counts)
dim(filt_counts)
sce <- fitGAM(counts = as.matrix(filt_counts), sds = curves)
plotGeneCount(curves, counts, clusters = clustering, models = sce)

# saveRDS(filt_counts, "rds/trajectory_filtCounts_withGOI.rds")
# saveRDS(sce, "rds/trajectory_fitGAM_withGOI.rds")
# saveRDS(curves, "rds/trajectory_curves_withGOI.rds")

filt_counts <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq_pub/slingshot/trajectory_filtCounts_withGOI.rds")
sce <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq_pub/slingshot/trajectory_fitGAM.rds")
curves <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq_pub/slingshot/trajectory_curves_withGOI.rds")

#define plot function
plot_differential_expression <- function(feature_id) {
  cowplot::plot_grid(plotGeneCount(curves, filt_counts, gene = feature_id[1], clusters = clustering, models = sce) + ggplot2::theme(legend.position = "none"), 
                     plotSmoothers(sce, as.matrix(counts), gene = feature_id[1]))
}

pseudotime_association <- associationTest(sce, lineages = TRUE)
# pseudotime_association$fdr <- p.adjust(pseudotime_association$pvalue, method = "fdr")
# pseudotime_association$fdr_1 <- p.adjust(pseudotime_association$pvalue_1, method = "fdr")
# pseudotime_association$fdr_2 <- p.adjust(pseudotime_association$pvalue_2, method = "fdr")
pseudotime_association <- pseudotime_association[order(pseudotime_association$pvalue), ]
pseudotime_association$feature_id <- rownames(pseudotime_association)

write.csv(pseudotime_association, "rds/association_test.csv")

feature_id <- "nbis-gene-2-utr"
plot_differential_expression(feature_id)


#PLOTTING 

cols <- c("0" = "#F8766D", "2" = "#D89000", "3" = "#C09B00", "5" = "#7CAE00", "6" = "#39B600", "10" = "#00BFC4", "11" = "#00BAE0", "12" = "#00B0F6", "13" = "#9590FF", "18" = "#FF6A98")


curves_plot <- plotGeneCount(curves, filt_counts, clusters = clustering, models = sce) + 
  ggplot2::theme(axis.text.x = element_text(size = 17)) + 
  ggplot2::theme(axis.text.y = element_text(size = 17)) + 
  #  labs(color = "WNV 5' UTR count") +
  #  theme(legend.position = c(0.2, 0.85)) +
  scale_x_continuous(limits = c(-12, 11)) +
  theme(axis.title = element_text(size = 17), legend.position = 'none') +
  scale_color_manual(values = cols)

exp_plot <- plotSmoothers(sce, as.matrix(counts), gene = feature_id[1], curvesCols = c("#7CAE00", "#00B0F6")) +
  ggplot2::theme(axis.text.x = element_text(size = 17)) + 
  ggplot2::theme(axis.text.y = element_text(size = 17)) +
  scale_color_manual(values = c("#7CAE00", "#00B0F6")) +
  theme(axis.title = element_text(size = 17))

plot <- curves_plot + exp_plot + patchwork::plot_layout(ncol = 2)



