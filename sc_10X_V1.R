library(DropletUtils)
library(singleCellTK)
library(Seurat)
library(tidyverse)
library(magrittr)
#library(usethis)

#generate sce input files from cell_ranger outs 
mg1_sce <- importCellRanger(sampleDirs = "/home/fitzmeyer/data_sets/cellranger_outs/wnv_mg1",
                            dataType = c("raw"),
                            matrixFileNames = "matrix.mtx.gz",
                            featuresFileNames = "features.tsv.gz",
                            barcodesFileNames = "barcodes.tsv.gz",
                            gzipped = TRUE)

#mg2_sce <- importCellRanger(sampleDirs = "/Users/emilyfitzmeyer/Documents/BC_WNV/scRNAseq/sc_outs/wnv_mg2",
#                            dataType = c("raw"),
#                            matrixFileNames = "matrix.mtx.gz",
#                           featuresFileNames = "features.tsv.gz",
#                            barcodesFileNames = "barcodes.tsv.gz",
#                            gzipped = TRUE)

#test emptyDrops performance
set.seed(100)
limit1 <- 40
#limit2 <- 40
out_1 <- emptyDrops(mg1_sce, lower = limit1, test.ambient = TRUE)
#out_2 <- emptyDrops(mg2_sce, lower = limit2, test.ambient = TRUE)

hist(out_1$PValue[out_1$Total <= limit1 & out_1$Total > 0],
     xlab="P-value", main="", col="palegreen") 

#hist(out_2$PValue[out_2$Total <= limit2 & out_2$Total > 0],
#     xlab="P-value", main="", col="palegreen") 

#subset sce object to retain only detected cells 
set.seed(100)
out_mg1 <- emptyDrops(counts(mg1_sce))
summary(out_mg1$FDR <= 0.001)
mg1_sce1 <- mg1_sce[,which(out_mg1$FDR <= 0.001)]

#set.seed(100)
#out_mg2 <- emptyDrops(counts(mg2_sce))
#summary(out_mg2$FDR <= 0.001)
#mg2_sce1 <- mg2_sce[,which(out_mg2$FDR <= 0.001)]

#if I define data="logcounts" in as.Seurat it throws a fit so I tried to define the log counts like so:
#counts <- assay(mg1_sce1, "counts")
#libsizes <- colSums(counts)
#size_factors <- libsizes/mean(libsizes)
#but this line throws a vector memory exhausted error that I wasn't able to fix by editing .Renviron:
#logcounts(mg1_sce1) <- log2(t(t(counts)/size_factors) + 1)
#mg1_seurat <- as.Seurat(mg1_sce1)

#reading data=NULL is a work around
#mg1_seurat <- as.Seurat(mg1_sce1, counts="counts", data=NULL)

mg1_seurat <- convertSCEToSeurat(mg1_sce1)
#mg2_seurat <- convertSCEToSeurat(mg2_sce1)

mg1_seurat[["percent.mt"]] <- PercentageFeatureSet(mg1_seurat, pattern = "^MT-")
#mg2_seurat[["percent.mt"]] <- PercentageFeatureSet(mg2_seurat, pattern = "^MT-")

VlnPlot(mg1_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#VlnPlot(mg2_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

mg1_plot1 <- FeatureScatter(mg1_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
mg1_plot2 <- FeatureScatter(mg1_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
mg1_plot2

#mg2_plot1 <- FeatureScatter(mg2_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
#mg2_plot2 <- FeatureScatter(mg2_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#mg2_plot2

mg1_seurat_VF <- FindVariableFeatures(mg1_seurat, selection.method = "vst", nfeatures = 2000)
mg1_top10 <- head(VariableFeatures(mg1_seurat_VF), 10)

#mg2_seurat_VF <- FindVariableFeatures(mg2_seurat, selection.method = "vst", nfeatures = 2000)
#mg2_top10 <- head(VariableFeatures(mg2_seurat_VF), 10)

#mg1 - Identification of highly variable features:
mg1_plot1_VF <- VariableFeaturePlot(mg1_seurat_VF, raster = FALSE)
mg1_plot2_VF <- LabelPoints(plot = mg1_plot1_VF, points = mg1_top10, repel = TRUE, xnudge = 0, ynudge = 0)
mg1_plot1_VF
mg1_plot2_VF

#mg2_plot1_VF <- VariableFeaturePlot(mg2_seurat_VF, raster = FALSE)
#mg2_plot2_VF <- LabelPoints(plot = mg2_plot1_VF, points = mg2_top10, repel = TRUE, xnudge = 0, ynudge = 0)
#mg2_plot1_VF
#mg2_plot2_VF

#clean variables (can avoid doing this by overwriting everything the way seurat tutorial does)
rm(mg1_plot1_VF, mg1_plot2_VF, mg1_plot2, mg1_plot1, mg1_seurat, mg1_sce1, mg1_sce, out_1, out_mg1)

#mg1 - scale data 
mg1_genes <- rownames(mg1_seurat_VF)
mg1_seurat_scaled <- ScaleData(mg1_seurat_VF, features = mg1_genes)

#run PCA on scaled data
mg1 <- RunPCA(mg1_seurat_scaled, features = VariableFeatures(object = mg1_seurat_scaled))
print(mg1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(mg1, dims = 1:2, reduction = "pca")
DimPlot(mg1, reduction = "pca")

#heatmaps
DimHeatmap(mg1, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(mg1, dims = 1:15, cells = 500, balanced = TRUE)

#determine dimensionality
mg1 <- JackStraw(mg1, num.replicate = 100)
mg1 <- ScoreJackStraw(mg1, dims = 1:20)
JackStrawPlot(mg1, dims = 1:15)
ElbowPlot(mg1)

#cluster cells
mg1 <- FindNeighbors(mg1, dims = 1:10)
mg1 <- FindClusters(mg1, resolution = 0.5)
head(Idents(mg1), 5)

#UMAP/tSNE
mg1 <- RunUMAP(mg1, dims = 1:10)
mg1 <- RunTSNE(mg1, dims = 1:10)
DimPlot(mg1, reduction = "umap")
DimPlot(mg1, reduction = "tsne")

#save object so it can be loaded back without having to rerun
saveRDS(mg1, file = "/home/fitzmeyer/data_sets/cellranger_outs/wnv_mg1_12dpi_seurat.rds")
