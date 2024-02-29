library(scLink)
library(tidyverse)

#read in genes (top 500 VFs) 
#and counts (transposed seurat object counts matrix) files 
#(generated in merge_Xdpi script):
#WORKING DIRECTORY FOR CCTSI-104

setwd("/home/fitzmeyer/rds")

genes <- readRDS("vf_4dpi_top500.rds")
count <- readRDS("4dpi_counts_matrix.rds")

count1 <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/scLink/4dpi_counts_matrix.rds")
ncol(count1)
nrow(count1)

count.norm <- sclink_norm(count, filter.genes = FALSE, gene.names = genes)

#alternatively you can just set filter.genes to 'TRUE' to filter for the top $n$ genes with largest average expression values
#count.norm = sclink_norm(count, scale.factor = 1e6, filter.genes = TRUE, n = 500)







# 
# 
# library(SingleCellExperiment)
# library(Seurat)
# library(scater)
# library(scran)
# library(igraph)
# library(pheatmap)
# library(mclust)
# library(tidyverse)
# 
# #InstallData("pbmc3k")
# pbmc <- LoadData("pbmc3k", type = "pbmc3k.final")
# view(pbmc@meta.data)
# 
# pbmc.sce <- as.SingleCellExperiment(pbmc)
# 
# #making a cell adjacency matrix using nearest neighbor/correlation:
# 
# graph_k10 <- scran::buildSNNGraph(pbmc.sce, k = 10, use.dimred = "PCA", type = "rank")
# 
# clust_k10_louvain <- igraph::cluster_louvain(graph_k10)$membership
# 
# table(clust_k10_louvain)
# 
# pbmc.sce$cluster_louvain_k10 <- factor(clust_k10_louvain)
# 
# scater::plotReducedDim(pbmc.sce, "UMAP", colour_by = "cluster_louvain_k10")
# 
# cols <- RColorBrewer::brewer.pal(n = 12, name = "Paired")
# 
# igraph::plot.igraph(
#   graph_k10, layout = layout_with_fr(graph_k10),
#   vertex.color = cols[clust_k10_louvain],
#   vertex.size = 5, vertex.label = NA, main = "Louvain"
# )






