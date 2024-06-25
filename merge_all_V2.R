library(singleCellTK)
library(tidyverse)
library(Seurat)
library(scCustomize)
library(gridExtra)
library(EnhancedVolcano)
library(clustree)
library(scater)
library(scran)
library(batchelor)
library(SeuratWrappers)

#12dpi
mg3.c.sce <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/mg3_c_QC.rds")

mg4.c.sce <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/mg4_c_QC.rds")

wnv.mg3.sce <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/wnv_mg3_QC.rds")

wnv.mg4.sce <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/wnv_mg4_QC.rds")

#4dpi
mg5.c.sce <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/mg5_c_QC.rds")

mg6.c.sce <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/mg6_c_QC.rds")

mg7.c.sce <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/mg7_c_QC.rds")

mg8.c.sce <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/mg8_c_QC.rds")

mg9.c.sce <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/mg9_c_QC.rds")

wnv.mg5.sce <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/wnv_mg5_QC.rds")

wnv.mg6.sce <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/wnv_mg6_QC.rds")

wnv.mg7.sce <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/wnv_mg7_QC.rds")

wnv.mg8.sce <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/wnv_mg8_QC.rds")

wnv.mg9.sce <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/wnv_mg9_QC.rds")

#subset:
#12dpi
mg3.c.sce <- subsetSCECols(mg3.c.sce, colData = c("decontX_contamination < 0.6",  
                                                  "scDblFinder_doublet_score < 0.9"))

mg4.c.sce <- subsetSCECols(mg4.c.sce, colData = c("decontX_contamination < 0.6",  
                                                  "scDblFinder_doublet_score < 0.9"))

wnv.mg3.sce <- subsetSCECols(wnv.mg3.sce, colData = c("decontX_contamination < 0.6",  
                                                      "scDblFinder_doublet_score < 0.9"))

wnv.mg4.sce <- subsetSCECols(wnv.mg4.sce, colData = c("decontX_contamination < 0.6",  
                                                      "scDblFinder_doublet_score < 0.9"))

#4dpi
mg5.c.sce <- subsetSCECols(mg5.c.sce, colData = c("decontX_contamination < 0.6",  
                                                  "scDblFinder_doublet_score < 0.9"))

mg6.c.sce <- subsetSCECols(mg6.c.sce, colData = c("decontX_contamination < 0.6",  
                                                  "scDblFinder_doublet_score < 0.9"))

mg7.c.sce <- subsetSCECols(mg7.c.sce, colData = c("decontX_contamination < 0.6",  
                                                  "scDblFinder_doublet_score < 0.9"))

mg8.c.sce <- subsetSCECols(mg8.c.sce, colData = c("decontX_contamination < 0.6",  
                                                  "scDblFinder_doublet_score < 0.9"))

mg9.c.sce <- subsetSCECols(mg9.c.sce, colData = c("decontX_contamination < 0.6",  
                                                  "scDblFinder_doublet_score < 0.9"))

wnv.mg5.sce <- subsetSCECols(wnv.mg5.sce, colData = c("decontX_contamination < 0.6",  
                                                      "scDblFinder_doublet_score < 0.9"))

wnv.mg6.sce <- subsetSCECols(wnv.mg6.sce, colData = c("decontX_contamination < 0.6",  
                                                      "scDblFinder_doublet_score < 0.9"))

wnv.mg7.sce <- subsetSCECols(wnv.mg7.sce, colData = c("decontX_contamination < 0.6",  
                                                      "scDblFinder_doublet_score < 0.9"))

wnv.mg8.sce <- subsetSCECols(wnv.mg8.sce, colData = c("decontX_contamination < 0.6",  
                                                      "scDblFinder_doublet_score < 0.9"))

wnv.mg9.sce <- subsetSCECols(wnv.mg9.sce, colData = c("decontX_contamination < 0.6",  
                                                      "scDblFinder_doublet_score < 0.9"))


#convert to Seurat objects and save srt.rds files (have to do conversion on command line after update, idk why):
#12dpi
mg3c <- convertSCEToSeurat(mg3.c.sce)
saveRDS(mg3c, "/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/srt_rds/mg3c.rds")
mg4c <- convertSCEToSeurat(mg4.c.sce)
saveRDS(mg4c, "/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/srt_rds/mg4c.rds")
wnvMg3 <- convertSCEToSeurat(wnv.mg3.sce)
saveRDS(wnvMg3, "/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/srt_rds/wnvMg3.rds")
wnvMg4 <- convertSCEToSeurat(wnv.mg4.sce)
saveRDS(wnvMg4, "/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/srt_rds/wnvMg4.rds")
#4dpi
mg5c <- convertSCEToSeurat(mg5.c.sce)
saveRDS(mg5c, "/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/srt_rds/mg5c.rds")
mg6c <- convertSCEToSeurat(mg6.c.sce)
saveRDS(mg6c, "/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/srt_rds/mg6c.rds")
mg7c <- convertSCEToSeurat(mg7.c.sce)
saveRDS(mg7c, "/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/srt_rds/mg7c.rds")
mg8c <- convertSCEToSeurat(mg8.c.sce)
saveRDS(mg8c, "/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/srt_rds/mg8c.rds")
mg9c <- convertSCEToSeurat(mg9.c.sce)
saveRDS(mg9c, "/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/srt_rds/mg9c.rds")
wnvMg5 <- convertSCEToSeurat(wnv.mg5.sce)
saveRDS(wnvMg5, "/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/srt_rds/wnvMg5.rds")
wnvMg6 <- convertSCEToSeurat(wnv.mg6.sce)
saveRDS(wnvMg6, "/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/srt_rds/wnvMg6.rds")
wnvMg7 <- convertSCEToSeurat(wnv.mg7.sce)
saveRDS(wnvMg7, "/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/srt_rds/wnvMg7.rds")
wnvMg8 <- convertSCEToSeurat(wnv.mg8.sce)
saveRDS(wnvMg8, "/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/srt_rds/wnvMg8.rds")
wnvMg9 <- convertSCEToSeurat(wnv.mg9.sce)
saveRDS(wnvMg9, "/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/srt_rds/wnvMg9.rds")

#START HERE!
#read srt RDS files in if restarting srt workflow
#12dpi
mg3c <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/srt_rds/mg3c.rds")
mg4c <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/srt_rds/mg4c.rds")
wnvMg3 <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/srt_rds/wnvMg3.rds")
wnvMg4 <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/srt_rds/wnvMg4.rds")
#4dpi
mg5c <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/srt_rds/mg5c.rds")
mg6c <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/srt_rds/mg6c.rds")
mg7c <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/srt_rds/mg7c.rds")
mg8c <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/srt_rds/mg8c.rds")
mg9c <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/srt_rds/mg9c.rds")
wnvMg5 <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/srt_rds/wnvMg5.rds")
wnvMg6 <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/srt_rds/wnvMg6.rds")
wnvMg7 <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/srt_rds/wnvMg7.rds")
wnvMg8 <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/srt_rds/wnvMg8.rds")
wnvMg9 <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/srt_rds/wnvMg9.rds")

#Seurat QC function
srt_QC <- function(x) {
  x[["percent_mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")
  x <- subset(x, subset = nFeature_RNA > 100 &
                nFeature_RNA < 2500 &
                percent_mt < 25)
}

#Seurat QC
#12dpi
mg3c <- srt_QC(mg3c)
mg4c <- srt_QC(mg4c)
wnvMg3 <- srt_QC(wnvMg3)
wnvMg4 <- srt_QC(wnvMg4)
#4dpi
mg5c <- srt_QC(mg5c)
mg6c <- srt_QC(mg6c)
#mg7c <- srt_QC(mg7c)
mg8c <- srt_QC(mg8c)
mg9c <- srt_QC(mg9c)
#wnvMg5 <- srt_QC(wnvMg5)
wnvMg6 <- srt_QC(wnvMg6)
wnvMg7 <- srt_QC(wnvMg7)
wnvMg8 <- srt_QC(wnvMg8)
wnvMg9 <- srt_QC(wnvMg9)

#Seurat normalization 
#12dpi
mg3c <- NormalizeData(mg3c)
mg4c <- NormalizeData(mg4c)
wnvMg3 <- NormalizeData(wnvMg3)
wnvMg4 <- NormalizeData(wnvMg4)
#4dpi
mg5c <- NormalizeData(mg5c)
mg6c <- NormalizeData(mg6c)
#mg7c <- NormalizeData(mg7c)
mg8c <- NormalizeData(mg8c)
mg9c <- NormalizeData(mg9c)
#wnvMg5 <- NormalizeData(wnvMg5)
wnvMg6 <- NormalizeData(wnvMg6)
wnvMg7 <- NormalizeData(wnvMg7)
wnvMg8 <- NormalizeData(wnvMg8)
wnvMg9 <- NormalizeData(wnvMg9)

#merge, retaining sample ID (add sample condition for later use):
merged.srt <- merge(mg3c, y = c(mg4c, mg5c, mg6c, mg8c, mg9c, wnvMg3, wnvMg4, wnvMg6, wnvMg7, wnvMg8, wnvMg9),
                    add.cell.ids = c("mg3c_mock_12dpi_R3", "mg4c_mock_12dpi_R3", "mg5c_mock_4dpi_R4", "mg6c_mock_4dpi_R4", "mg8c_mock_4dpi_R6", 
                                     "mg9c_mock_4dpi_R6", "wnvMg3_wnv_12dpi_R3", "wnvMg4_wnv_12dpi_R3", "wnvMg6_wnv_4dpi_R5", "wnvMg7_wnv_4dpi_R5", 
                                     "wnvMg8_wnv_4dpi_R5", "wnvMg9_wnv_4dpi_R5"),
                    project = "mg_all", merge.data = TRUE)

#view(merged.srt@meta.data)
merged.srt$sample <- rownames(merged.srt@meta.data)
merged.srt@meta.data <- separate(merged.srt@meta.data, col = 'sample', into = c("sample", "condition", "dpi", "run", "barcode"),
                                 sep = "_")

#saveRDS(merged.srt, file = "/Users/emilyfitzmeyer/Desktop/alldpi_nmn_oldNorm_mergedObj.rds")
merged.srt <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq_pub/mergeDPI_oldNorm/alldpi_nmn_oldNorm_mergedObj.rds")

#Seurat QC:
# merged.srt[["percent_mt"]] <- PercentageFeatureSet(merged.srt, pattern = "^MT-")
# merged.srt <- subset(merged.srt, subset = nFeature_RNA > 100 & 
#                        nFeature_RNA < 2500 & 
#                        percent_mt < 25)
# 
# merged.srt <- srt_QC(merged.srt)

#generate and transpose counts matrix for scLink:
# split.srt1 <- SplitObject(merged.srt, split.by = "condition")
# wnv.srt1 <- split.srt1$wnv
# unique(wnv.srt1$sample)
# counts.matrix <- wnv.srt1@assays$RNA@counts
# t.4dpi.counts.matrix <- t(counts.matrix)
# saveRDS(t.4dpi.counts.matrix, "/Users/emilyfitzmeyer/Desktop/4dpi_wnv_counts_matrix.rds")

# wnv.srt1 <- NormalizeData(wnv.srt1)
# wnv.srt1 <- FindVariableFeatures(wnv.srt1)

#making gene of interest list for scLink:
# vf.wnv4dpi.srt <- VariableFeatures(wnv.srt1)
# vf.4dpi.top500 <- head(vf.wnv4dpi.srt, 500)
# saveRDS(vf.4dpi.top500, "/Users/emilyfitzmeyer/Desktop/vf_wnv4dpi_top500.rds")

#norm:
merged.srt <- NormalizeData(merged.srt)
merged.srt <- FindVariableFeatures(merged.srt)

#TEMP BEC STUFF
# merged.srt <- SplitObject(merged.srt, split.by = "run")
# merged.srt <- map(merged.srt, as.SingleCellExperiment)
# 
# r3 <- merged.srt$R3
# r4 <- merged.srt$R4
# r5 <- merged.srt$R5
# r6 <- merged.srt$R6
# 
# r3.dec <- modelGeneVar(r3)
# r4.dec <- modelGeneVar(r4)
# r5.dec <- modelGeneVar(r5)
# r6.dec <- modelGeneVar(r6)
# 
# combined.dec <- combineVar(r3.dec, r4.dec, r5.dec, r6.dec)
# chosen.hvgs <- getTopHVGs(combined.dec, n=5000)
# 
# combined <- correctExperiments(A=r3, B=r4, C=r5, D=r6, PARAM=NoCorrectParam())
# combined <- runPCA(combined, subset_row=chosen.hvgs)
# combined <- runUMAP(combined, dimred="PCA")
# plotUMAP(combined, colour_by="run")
# 
# mnn.out <- fastMNN(combined, batch=combined$run, subset.row=chosen.hvgs)
# str(reducedDim(mnn.out, "corrected"))
# 
# mnn.out <- runUMAP(mnn.out, dimred="corrected")
# plotUMAP(mnn.out, colour_by="batch")

#BEC w/ seurat wrappers:
# DefaultAssay(merged.srt) <- "RNA"
# merged.srt <- RunFastMNN(object.list = SplitObject(merged.srt, split.by = "run"))

#scale:
#srt.genes <- rownames(merged.srt)
#optional "features = srt.genes" argument
#optional vars.to.regress = "percent.mt" argument
merged.srt <- ScaleData(merged.srt)


#PCA
merged.srt <- RunPCA(merged.srt)
#elbow plot (choose PCs)
#ElbowPlot(merged.srt, ndims = 40)

#find neighbors:
merged.srt <- FindNeighbors(merged.srt, reduction = "pca", dims = 1:35)

#test clustering resolution
#merged.srt <- FindClusters(merged.srt, resolution = c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6))
#clustree(merged.srt)
#setwd("/Users/emilyfitzmeyer/Desktop/")
#ggsave("clustree_alldpi.png", plot = last_plot(), device = png(), scale = 1, width = 12.5, height = 10.5, dpi = 300)
#dev.off()

#find clusters:
merged.srt <- FindClusters(merged.srt, resolution = 0.6)
merged.srt <- RunUMAP(merged.srt, reduction = "pca", dims = 1:35)
#save merged object:
saveRDS(merged.srt, file = "/Users/emilyfitzmeyer/Desktop/alldpi_oldNorm_nmn_res0.6.rds")

#Read in merged.srt
merged.srt <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq_pub/mergeDPI_oldNorm/alldpi_oldNorm_nmn_res0.6.rds")
#merged.srt <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq_pub/mt_regress_test/alldpi_oldNorm_nmn_mtRegress_res0.6.rds")

#for terminal
#merged.srt <- readRDS("/home/fitzmeyer/rds/4dpi_res0.6.rds")

#rename clusters
new_idents <- c("EC-like-1", "1", "EC-like-2", "ISC/EB", "4", "EC", "EC-like-3", "VM-1", "8", "VM-2", "ISC/EB-prol", "cardia-1", "EE", "cardia-2", "HC-1", "MT", "HC-2", "17", "cardia-prol", "19")
merged.srt <- Rename_Clusters(merged.srt, new_idents = new_idents)

#add cell type info to metadata
ident_num <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19")
cell_type <- c("EC-like-1", "1", "EC-like-2", "ISC/EB", "4", "EC", "EC-like-3", "VM-1", "8", "VM-2", "ISC/EB-prol", "cardia-1", "EE", "cardia-2", "HC-1", "MT", "HC-2", "17", "cardia-prol", "19")
tib <- tibble(ident_num, cell_type)
merged.srt <- Add_Sample_Meta(merged.srt, meta_data = tib, join_by_seurat = "seurat_clusters", join_by_meta = "ident_num")

#FIND CLUSTER MARKERS
#===============================================================================
#set default assay to "RNA" so DE analyses use raw counts:
DefaultAssay(merged.srt) <- "RNA"
split.srt <- SplitObject(merged.srt, split.by = "dpi")
srt_4dpi<- split.srt$"4dpi"
srt_12dpi <- split.srt$"12dpi"
srt_4dpi <- SplitObject(srt_4dpi, split.by = "condition")
srt_12dpi <- SplitObject(srt_12dpi, split.by = "condition")

#Find all cluster markers for wnv+ samples:
markers_all <- FindAllMarkers(srt_12dpi$wnv, 
                              logfc.threshold = 0.25,
                              test.use = "wilcox")

write.csv(markers_all, file = "/Users/emilyfitzmeyer/Desktop/all_markers_wnv12dpi_wilcox.csv")

#and mock samples:
markers_mock <- FindAllMarkers(srt_12dpi$wnv, 
                               logfc.threshold = 0.25,
                               test.use = "wilcox")

write.csv(markers_mock, file = "/Users/emilyfitzmeyer/Desktop/wnv12_all_markers_wilcox.csv")

#identify all/conserved cell type markers:
cnsrv.markers <- FindConservedMarkers(merged.srt, 
                                      ident.1 = 5,
                                      grouping.var = "condition",
                                      slot = "data",
                                      test.use = "wilcox",
                                      verbose = FALSE)

write.csv(cnsrv.markers, file = "/Users/emilyfitzmeyer/Desktop/cl5_alldpi_conservedMarkers_mtRegress_wilcox.csv")

#find markers that differentiate clusters 
cl_markers <- FindMarkers(srt_4dpi$wnv,
                          ident.1 = 14,
                          slot = "data",
                          test.use = "wilcox",
                          verbose = FALSE)

write.csv(cl_markers, file = "/Users/emilyfitzmeyer/Desktop/cl14_v_all_markers_mockOnly_wilcox.csv")

#DESeq2 pseudo-bulk 
#===============================================================================
#srt.obj <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/integrated_mg3mg4.rds")
#srt.obj <- readRDS("/Users/emilyfitzmeyer/Desktop/temp.rds")

#create a column that contains both sample and condition information:
srt.obj$samp_cond <- paste(srt.obj$condition, srt.obj$sample, sep = "_")

#aggregate counts to the sample level:
counts <- AggregateExpression(srt.obj, group.by = "samp_cond",
                              assays = 'RNA',
                              slot = "counts",
                              return.seurat = FALSE)

#get counts matrix for DESeq2:
counts <- counts$RNA

#build DESeq2 dataset:
#generate sample level metadata:
colData <- data.frame(samples = colnames(counts), row.names = NULL)
colData <- colData %>%
  mutate(condition = ifelse(grepl('wnv', samples), 'wnv', 'mock')) %>%
  column_to_rownames(var = 'samples')

#create DESeq2 object:
dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, 
                                      colData = colData, 
                                      design = ~ condition)

#filter:
dds.keep <- rowSums(counts(dds)) >=10
dds <- dds[dds.keep,]

#run DESeq2:
dds <- DESeq2::DESeq(dds)

#generate results object:
DESeq2::resultsNames(dds)
dds.results <- DESeq2::results(dds, name = "condition_wnv_vs_mock")

#save results in various forms:
#saveRDS(dds.results, file = "/Users/emilyfitzmeyer/Desktop/4dpi_DESeq2_results_0.6.rds", compress = TRUE)
write.csv(dds.results, file = "/Users/emilyfitzmeyer/Desktop/temp.csv")

dds.results <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/merged_rds/4dpi_DESeq2_results_0.6.rds")
results <- read.csv("/Users/emilyfitzmeyer/Desktop/scRNAseq/plots/4dpi_labChange_volcanoPlotValues.csv")

#cutoff varies based on results file
EnhancedVolcano(results, 
                lab = results$gene_ID,
                title = "WNV vs. Mock DEGs 4dpi",
                subtitle = NULL,
                caption = NULL,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.000035,
                drawConnectors = TRUE,
                arrowheads = FALSE,
                labSize = 5,
                axisLabSize = 23,
                legendIconSize = 3)

setwd("/Users/emilyfitzmeyer/Desktop/")
ggsave("PB_DEGs_4dpi.png", plot = last_plot(), device = png(), scale = 1, width = 6, height = 6, dpi = 300)
dev.off()


#VISUALIZATION
#===============================================================================

#split obj for feature plots and marker identification
split.srt <- SplitObject(merged.srt, split.by = "dpi")
srt_4dpi <- split.srt$"4dpi"
srt_12dpi <- split.srt$"12dpi"
srt_4dpi <- SplitObject(srt_4dpi, split.by = "condition")
srt_12dpi <- SplitObject(srt_12dpi, split.by = "condition")


#VISUALIZE CLUSTERS:
#===============================================================================
#Examine clusters in each condition
clusters_alldpi <- DimPlot(merged.srt, reduction = "umap", label = TRUE, repel = TRUE, split.by = "dpi", pt.size = 0.2) +
  theme(axis.text.x = element_text(size = 15)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(legend.position = 'none')

VlnPlot(merged.srt, features = "nbis-gene-2-utr")
FeaturePlot(merged.srt, features = "gene10632")

ggsave("clusters_alldpi.png", plot = last_plot(), device = png(), scale = 1, width = 8.5, height = 5, dpi = 300)
dev.off()

#Examine clusters in each condition - NO SPLIT
clusters_4dpi <- DimPlot(merged.srt, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.2) +
  theme(axis.text.x = element_text(size = 15)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(legend.position = 'none')

ggsave("clusters_4dpi.png", plot = last_plot(), device = png(), scale = 1, width = 7, height = 6, dpi = 300)
dev.off()

#show sample and condition induced impact on cell grouping:
#===============================================================================
#group by sample or run
sample_plot <- DimPlot(srt_12dpi, reduction = 'umap', group.by = 'sample') &
  theme(axis.text.x = element_text(size = 17)) &
  theme(axis.text.y = element_text(size = 17)) &
  ggtitle("sample grouping")
ggsave("run_grouping.png", plot = last_plot(), device = png(), scale = 1, width = 6.5, height = 5, dpi = 300)
dev.off()
#group by condition or dpi
condition_plot <- DimPlot(merged.srt, reduction = 'umap', group.by = 'condition',
                          cols = c('deepskyblue', 'salmon')) &
  theme(axis.text.x = element_text(size = 17)) &
  theme(axis.text.y = element_text(size = 17)) &
  ggtitle("condition grouping")
ggsave("condition_grouping.png", plot = last_plot(), device = png(), scale = 1, width = 6.5, height = 5, dpi = 300)
dev.off()

#subset clusters if needed:
#===============================================================================
cl4.srt <- subset(merged.srt, idents = 4)
DimPlot(cl4.srt, reduction = 'umap', group.by = 'sample')
VlnPlot(cl4.srt, group.by = 'sample', features = "gene4107")
VlnPlot(merged.srt, features = "gene4107")
FeaturePlot(merged.srt, features = "gene14035")

#WNV vRNA load calculation and visualization
#===============================================================================
mock.percent.whole <- Percent_Expressing(split.srt$mock, features = "nbis-gene-2-utr", entire_object = TRUE)
wnv.percent.whole <- Percent_Expressing(split.srt$wnv, features = "nbis-gene-2-utr", entire_object = TRUE)
wnv.percent.cluster <- Percent_Expressing(split.srt$wnv, features = "nbis-gene-2-utr")
wnv.expLevel.whole <- AverageExpression(split.srt$wnv, features = "nbis-gene-2-utr", group.by = "orig.ident")

#plot presence of WNV 5' UTR:
#===============================================================================
WNV.utr.plot <- FeaturePlot(merged.srt, features = "nbis-gene-2-utr", split.by = "condition", pt.size = 0.5) &
  theme(axis.text.x = element_text(size = 17)) &
  theme(axis.text.y = element_text(size = 17)) &
  theme(legend.position = c(0.9, 0.85))


#WNV vRNA level VLN plots:
setwd("/Users/emilyfitzmeyer/Desktop/")

my_levels <- c("VM-2", "VM-1", "MT", "ISC/EB-prol", "ISC/EB", "HC-2", "HC-1", "EE", "EC-like-3", "EC-like-2", "EC-like-1", "EC", "cardia-2", "cardia-1", "cardia-prol", "19", "17", "8", "4", "1")
factor(Idents(srt_4dpi$wnv), levels = my_levels)
Idents(srt_4dpi$wnv) <- factor(Idents(srt_4dpi$wnv), levels = my_levels)

vrna_level_4dpi_vlnPlot <- VlnPlot(srt_4dpi$wnv, features = "nbis-gene-2-utr", idents = my_levels, cols = c("#CC99FF","#CC99FF","#CC99FF","#CC99FF","#CC99FF","#CC99FF","#CC99FF","#CC99FF","#CC99FF","#CC99FF","#CC99FF","#CC99FF","#CC99FF","#CC99FF","#CC99FF","#CC99FF","#CC99FF","#CC99FF","#CC99FF", "#CC99FF")) +
  theme(axis.text.x = element_text(size = 15)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(legend.position = 'none') +
  #scale_x_discrete(labels = NULL) +
  labs(x = NULL)

write.csv(vrna_level_4dpi_vlnPlot$data, file = "/Users/emilyfitzmeyer/Desktop/vRNA_vlnplot_4dpi.csv")

factor(Idents(srt_12dpi$wnv), levels = my_levels)
Idents(srt_12dpi$wnv) <- factor(Idents(srt_12dpi$wnv), levels = my_levels)

vrna_level_12dpi_vlnPlot <- VlnPlot(srt_12dpi$wnv, features = "nbis-gene-2-utr", idents = my_levels, cols = c("darkslategray3","darkslategray3","darkslategray3","darkslategray3","darkslategray3","darkslategray3","darkslategray3","darkslategray3","darkslategray3","darkslategray3","darkslategray3","darkslategray3","darkslategray3","darkslategray3","darkslategray3","darkslategray3","darkslategray3")) + 
  theme(axis.text.x = element_text(size = 15)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(legend.position = 'none') +
  #scale_x_discrete(labels = NULL) +
  labs(x = NULL)

write.csv(vrna_level_12dpi_vlnPlot$data, file = "/Users/emilyfitzmeyer/Desktop/vRNA_vlnplot_12dpi.csv")

ggsave("12dpi_wnv_pct_cluster.png", plot = last_plot(), device = png(), scale = 1, width = 5.5, height = 3.5, dpi = 300)
dev.off()

#get plot data:
write.csv(plot_temp$data, file = "/Users/emilyfitzmeyer/Desktop/knwn_CLs_WNV_vRNA.csv")

#CLUSTER MARKER DOTPLOT
#===============================================================================
dotplot_genes <- c("gene11957", "gene316", "gene769", "gene10632", "gene2985", "gene9804", "gene11104", "gene2439", 
                   "gene452", "gene7512", "gene1416", "gene13104", "gene2167", "gene12446", "gene9049", "gene3360", "gene3283", "gene11820", "gene6403", 
                   "gene6861", "gene13589", "gene14180", "gene5110", "gene4641", "gene12764",
                   "gene10257", "gene927", "gene7675")

dotplot_break_labels <- c("POU2F1", "PLA2G6", "AGBL5", "PROX1", "ACTB", "Mlc2", "Mhc", "NIMB2", "SPARC", "pebIII_CPIJ002629", "pebIII_CPIJ002609", 
                          "klumpfuss", "PCNA", "C-type lysozyme", "Sugar_tr_CPIJ019592", "Sugar_tr_CPIJ014327", "Sugar_tr_CPIJ011910",
                          "Sugar_tr_CPIJ012675", "Sugar_tr_CPIJ012678", "chitin-binding_CPIJ004734", "chitin-binding_CPIJ004734",
                          "serine_protease_CPIJ006568", "serine_protease_CPIJ015103",
                          "serine_protease_CPIJ007079", "serine_protease_KDR18614", 
                          "Mal-B2", "Mal-A4", "irk-2")

#arrange "idents" alphabetically for plot (this method puts the numbered clusters in a weird order)
# levels(merged.srt@meta.data$RNA_snn_res.0.6)
# 
# Idents(merged.srt) <- "cell_type"
# Idents(merged.srt) <- factor(merged.srt@active.ident, sort(levels(merged.srt@active.ident)))

#reorder y axis to list cell types alphabetically/group types together (no number weirdness)
#For some reason this approach stopped working with allDPI merge 
my_levels <- c("VM-2", "VM-1", "MT", "ISC/EB-prol", "ISC/EB", "HC-2", "HC-1", "EE", "EC-like-3", "EC-like-2", "EC-like-1", "EC", "cardia-2", "cardia-1", "cardia-prol", "19", "17", "8", "4", "1")
factor(Idents(merged.srt), levels = my_levels)
Idents(merged.srt) <- factor(Idents(merged.srt), levels = my_levels)

#back to old approach for ordering clusters:
Idents(merged.srt) <- "cell_type"
Idents(merged.srt) <- factor(merged.srt@active.ident, sort(levels(merged.srt@active.ident)))

 DotPlot_scCustom(merged.srt, features = dotplot_genes) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = dotplot_break_labels)

setwd("/Users/emilyfitzmeyer/Desktop/")
ggsave("cl_marker_dotplot_alldpi.png", plot = last_plot(), device = png(), scale = 1, width = 9, height = 6.5, dpi = 300)
dev.off()

#view cell cycle features
#===============================================================================
#S phase
FeaturePlot(merged.srt, features = c("gene7803", "gene1526", "gene2166", "gene577", "gene9562"))
#G2/M phase
FeaturePlot(merged.srt, features = c("gene13155", "gene9499", "gene7697", "gene4233", "gene11073", "gene13413", "gene13452", "gene5102", "gene9506", "gene6438", "gene7563", "gene9947"))


#mito gene enriched clusters dotplot
#===============================================================================
dotplot_genes <- c("nbisL1-trna-9", "MT-nbis-gene-3", "MT-nbis-gene-2", "nbisL1-trna-7", "nbisL1-trna-10",
                   "nbisL1-trna-6", "nbisL1-trna-16", "nbisL1-trna-3", "nbisL1-trna-17", "MT-nbis-gene-5",
                   "nbisL1-trna-20", "nbisL1-trna-19", "nbis-gene-2-utr")

dotplot_break_labels <- c("MT-tRNA-Asp", "COX3", "ATP6", "MT-tRNA-Leu", "MT-tRNA-Gly",
                          "MT-tRNA-Tyr", "MT-tRNA-Phe", "MT-tRNA-Met", "MT-tRNA-His", "CYTB",
                          "MT-tRNA-Ser", "MT-tRNA-Pro", "WNV 5' UTR")

#reorder y axis to list cell types alphabetically/group types together (no number weirdness)
#my_levels <- c("VM-2", "VM-1", "ISC/EB-prol", "ISC/EB", "HC-2", "HC-1", "EE", "EC-like-2", "EC-like-1", "EC", "cardia-2", "cardia-1", "15", "5", "3", "1")
#fix outdated levels 
factor(Idents(split.srt$wnv), levels = my_levels)
Idents(split.srt$wnv) <- factor(Idents(split.srt$wnv), levels = my_levels)

DotPlot_scCustom(split.srt$wnv, features = dotplot_genes) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = dotplot_break_labels)

ggsave("cl3_marker_dotplot.png", plot = last_plot(), device = png(), scale = 1, width = 6, height = 4, dpi = 300)
dev.off()



#ISC/EB clusters dotplot - Irrelevant with dpi merge
#===============================================================================
dotplot_genes <- c("gene6060", "gene4986", "gene10089",	"gene6115",	"gene6292",	"gene1252",	"gene6564",
                   "gene2793", "gene6552", "gene5534", "gene4419")

dotplot_break_labels <- c("DUF4803", "gene4986", "cecropin", "cecropin", "CUTA", "GSTD11", "awd",
                          "NHP2L1", "histone H2A", "RAN", "BIRC5")

DotPlot_scCustom(merged.srt, features = dotplot_genes) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = dotplot_break_labels)

ggsave("ISC_marker_dotplot_4dpi.png", plot = last_plot(), device = png(), scale = 1, width = 6, height = 4, dpi = 300)
dev.off()


#VISUALIZE SPECIFIC GENES (correlation, expression, UMAP expression)
#===============================================================================
FeatureScatter(split.srt$wnv, feature1 = "nbis-gene-2-utr", feature2 = "gene9710", plot.cor = FALSE) + 
  theme(legend.position = "none") +
  labs(y = "ninjurin", x = "WNV 5' UTR")

VlnPlot(merged.srt, features = "gene2167", group.by = "condition") +
  theme(legend.position = "none") +
  ggtitle("PCNA")

#ISC plots
FeaturePlot(merged.srt, features = "gene13104", pt.size = 0.5) +
  theme(axis.text.x = element_text(size = 17)) +
  theme(axis.text.y = element_text(size = 17)) +
  scale_x_continuous(limits = c(-5, 4)) +
  scale_y_continuous(limits = c(-11, -6)) +
  theme(legend.position = c(0.9, 0.23)) +
  theme(legend.text = element_text(size = 17)) +
  ggtitle("Klu")
ggsave("klu_4dpi.png", plot = last_plot(), device = png(), scale = 1, width = 6, height = 5, dpi = 300)
dev.off()

#PROX1 plots
FeaturePlot(merged.srt, features = "gene10632") +
  theme(axis.text.x = element_text(size = 17)) +
  theme(axis.text.y = element_text(size = 17)) +
  ggtitle("PROX1")
ggsave("PROX1_4dpi.png", plot = last_plot(), device = png(), scale = 1, width = 6, height = 6, dpi = 300)
dev.off()
VlnPlot(merged.srt, features = "gene10632") + 
  theme(legend.position = 'none') + 
  theme(axis.text.x = element_text(size = 17)) +
  theme(axis.text.y = element_text(size = 17)) +
  ggtitle("PROX1")
ggsave("PROX1_vlnPlot_4dpi.png", plot = last_plot(), device = png(), scale = 1, width = 6, height = 6, dpi = 300)
dev.off()

#percent MT by cluster
#===============================================================================
temp <- group_by(merged.srt@meta.data, seurat_clusters)
list_temp <- group_split(temp)
mt_mean <- function(x) {
  mt_mean <- mean(x[["percent_mt"]])
}
lapply(list_temp, mt_mean)


#CLUSTER PROPORTIONS
#===============================================================================
cluster_proportions <- prop.table(table(Idents(merged.srt)))
as.data.frame(cluster_proportions)


#identifying pctExp>75 genes for COG analysis
#===============================================================================
#my adaptation - chatGPT version wasn't quite doing what I wanted
#note - my version needs to be run on the server bc it's chonky (how chonky? Run it in a screen.)
clusters <- unique(merged.srt$seurat_clusters)
gene_matrices <- list()

for (cluster in clusters) {
  subset_seurat <- subset(merged.srt, ident = cluster)
  cluster_genes <- rownames(subset_seurat@assays$RNA@data)
  pct_exp <- Percent_Expressing(subset_seurat, features = cluster_genes, entire_object = TRUE)
  top_75pct_genes <- subset(pct_exp, pct_exp$All_Cells >= 75)
  gene_matrices[[cluster]] <- top_75pct_genes
  print("next cluster")
}

saveRDS(gene_matrices, file = "/home/fitzmeyer/rds/gene_matrices_4dpi.rds")
#fetch file from server and proceed on local:
gene_matrices <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq_pub/COG/gene_matrices_4dpi.rds")

cl_nums <- names(gene_matrices)
cl_nums <- paste0(cl_nums, "pctExp75.csv")
setwd("/Users/emilyfitzmeyer/Desktop/")
#add colname to list elements
colnames <- c("gene_ID", "percent_expressing")
gene_matrices <- lapply(gene_matrices, rownames_to_column)
for (i in seq_along(gene_matrices)){
  colnames(gene_matrices[[i]]) <- colnames
}
#write out .csv files
for(i in 1:length(gene_matrices)){
  write.csv(gene_matrices[[i]], file = cl_nums[i], row.names = FALSE)
}



#percentExpression key immune genes:
#===============================================================================
immn.gene.pctExp <- Percent_Expressing(merged.srt, assay = 'RNA', split_by = "condition", features = c("gene1501", "gene7624", "gene10322", "gene11340", "gene3044",
                                                                                                       "gene6258", "gene8444", "gene11364", "gene6091", "gene6092",
                                                                                                       "gene5869", "gene5870", "gene10292", "gene10293", "gene3397",
                                                                                                       "gene14111", "gene2893", "gene3913", "gene2450", "gene5218", 
                                                                                                       "gene13590", "gene3988", "gene3989", "gene4535", "gene14091"),
                                       entire_object = TRUE)

#averageExpression key immune genes:
#===============================================================================

immn.gene.avgExp <- AverageExpression(merged.srt, assays = 'RNA', features = c("gene1501", "gene7624", "gene10322", "gene11340", "gene3044",
                                                                               "gene6258", "gene8444", "gene11364", "gene6091", "gene6092",
                                                                               "gene5869", "gene5870", "gene10292", "gene10293", "gene3397",
                                                                               "gene14111", "gene2893", "gene3913", "gene2450", "gene5218", 
                                                                               "gene13590", "gene3988", "gene3989", "gene4535", "gene14091"), 
                                      group.by = "condition")



#FIGURE 6F-G
split.srt <- SplitObject(merged.srt, split.by = "condition")

scatter_plot <- FeatureScatter(split.srt$wnv, feature1 = "nbis-gene-2-utr", feature2 = "gene11340", plot.cor = FALSE, group.by = "orig.ident", cols = "blue") + 
  theme(legend.position = "none") +
  labs(y = "IMD", x = "WNV 5' UTR") +
  theme(axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18)) +
  theme(axis.text.x = element_text(size = 17), axis.text.y = element_text(size = 17))

ggsave("alldpi_DOME_v_WNVsc.png", plot = scatter_plot, device = png(), scale = 1, width = 4, height = 5, dpi = 300)
dev.off()

vln_plot <- VlnPlot(merged.srt, features = "gene11340", group.by = "condition", pt.size = 0.5) +
  theme(legend.position = "none") +
  ggtitle("IMD") +
  theme(axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18)) +
  theme(axis.text.x = element_text(size = 17), axis.text.y = element_text(size = 17)) +
  scale_x_discrete(labels = c("Mock", "WNV")) +
  labs(x = NULL)

ggsave("alldpi_DOME_v_WNVvln.png", plot = vln_plot, device = png(), scale = 1, width = 2, height = 5, dpi = 300)
dev.off()


FeaturePlot(merged.srt, features = "nbis-gene-2-utr")
