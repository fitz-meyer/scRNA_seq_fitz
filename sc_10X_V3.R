library(Seurat)
library(singleCellTK)
library(tidyverse)
library(clustree)

#srt__________________________________________________________________________________________________
sce <- readRDS("/home/fitzmeyer/data_sets/scRNAseq_rds/aggr_mg3mg4_QC.rds")

plotDecontXResults(sce, reducedDimName = "decontX_UMAP")

#subset based on QC scores
sce <- subsetSCECols(sce, colData = c("decontX_contamination < 0.6",  
                                      "scDblFinder_doublet_score < 0.9"))

srt <- convertSCEToSeurat(sce)

srt[["percent_mt"]] <- PercentageFeatureSet(srt, pattern = "^MT-")
VlnPlot(srt, features = c("nFeature_RNA", "nCount_RNA", "percent_mt"), ncol = 3)

#plot1 <- FeatureScatter(srt, feature1 = "nCount_RNA", feature2 = "percent_mt")
#plot2 <- FeatureScatter(srt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

#subset based on feature values
srt <- subset(srt, subset = nFeature_RNA > 100 & nFeature_RNA < 2500 & percent_mt < 30)

#NORMALIZE
srt <- NormalizeData(srt)

srt <- FindVariableFeatures(srt, selection.method = "vst", nfeatures = 2000)
#srt_top10 <- head(VariableFeatures(srt), 10)

#srt_plot1_VF <- VariableFeaturePlot(srt, raster = FALSE)
#srt_plot2_VF <- LabelPoints(plot = srt_plot1_VF, points = srt_top10, repel = TRUE, 
#                            max.overlaps = 20, xnudge = 0, ynudge = 0)

#SCALE
#srt_genes <- rownames(srt)
#optional "features = srt_genes" argument
srt <- ScaleData(srt)
#saveRDS(srt, file = "/home/fitzmeyer/data_sets/cellranger_outs/srt_scaled_seurat.rds", compress = TRUE)

#run PCA
srt <- RunPCA(srt, features = VariableFeatures(object = srt))
#print(srt[["pca"]], dims = 1:5, nfeatures = 5)
#VizDimLoadings(srt, dims = 1:2, reduction = "pca")
#DimPlot(srt, reduction = "pca")

#ElbowPlot(srt)

view(srt@meta.data)

#change dims based on Elbow plot PCs 
srt <- FindNeighbors(srt, dims = 1:18)
srt <- FindClusters(srt, resolution = 0.7)
#head(Idents(srt), 5)

#clustree(srt)

srt <- RunUMAP(srt, dims = 1:18)
srt <- RunTSNE(srt, dims = 1:18)
DimPlot(srt, reduction = "umap")
DimPlot(srt, reduction = "tsne")

saveRDS(srt, file = "/home/fitzmeyer/data_sets/scRNAseq_rds/wnv_mg4_dimreduc.rds", compress = TRUE)

#srt <- readRDS("/home/fitzmeyer/data_sets/scRNAseq_rds/mg1_c_dimreduc.rds")

srt_markers <- FindAllMarkers(srt, only.pos = TRUE, 
                              min.pct = 0.25, logfc.threshold = 0.25, test.use = "DESeq2")

write.csv(srt_markers, file = "/home/fitzmeyer/data_sets/scRNAseq_rds/wnv_mg4_all_markers_DESeq2.csv")

FeaturePlot(srt, features = c("nbis-gene-2-utr"))

#mito gene feature plot:
#FeaturePlot(srt, features = c("nbis-gene-5", "nbisL1-trna-20", "nbisL1-trna-6", 
#                              "nbisL1-trna-3", "nbisL1-trna-7", "nbisL1-trna-9", 
#                              "nbis-gene-2", "nbis-gene-3", "nbisL1-trna-10", 
#                              "nbisL1-trna-16", "nbisL1-trna-17", "nbis-gene-4", 
#                              "nbisL1-trna-19")) 











