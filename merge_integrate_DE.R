library(singleCellTK)
library(tidyverse)
library(Seurat)
library(gridExtra)
library(EnhancedVolcano)

#read in SCE objects:
mg3_c <- importCellRanger(sampleDirs = "/Users/emilyfitzmeyer/Desktop/scRNAseq/cellranger_outs/mg3_c/",
                        dataType = c("filtered"),
                        matrixFileNames = "matrix.mtx.gz",
                        featuresFileNames = "features.tsv.gz",
                        barcodesFileNames = "barcodes.tsv.gz",
                        gzipped = TRUE)

mg4_c <- importCellRanger(sampleDirs = "/Users/emilyfitzmeyer/Desktop/scRNAseq/cellranger_outs/mg4_c/",
                        dataType = c("filtered"),
                        matrixFileNames = "matrix.mtx.gz",
                        featuresFileNames = "features.tsv.gz",
                        barcodesFileNames = "barcodes.tsv.gz",
                        gzipped = TRUE)

wnv_mg3 <- importCellRanger(sampleDirs = "/Users/emilyfitzmeyer/Desktop/scRNAseq/cellranger_outs/wnv_mg3/",
                        dataType = c("filtered"),
                        matrixFileNames = "matrix.mtx.gz",
                        featuresFileNames = "features.tsv.gz",
                        barcodesFileNames = "barcodes.tsv.gz",
                        gzipped = TRUE)

wnv_mg4 <- importCellRanger(sampleDirs = "/Users/emilyfitzmeyer/Desktop/scRNAseq/cellranger_outs/wnv_mg4/",
                        dataType = c("filtered"),
                        matrixFileNames = "matrix.mtx.gz",
                        featuresFileNames = "features.tsv.gz",
                        barcodesFileNames = "barcodes.tsv.gz",
                        gzipped = TRUE)

#estimate ambient RNA and identify doublets:
mg3_c <- runCellQC(mg3_c, sample = NULL,
                 algorithms = c("scDblFinder", "decontX"),
                 geneSetListLocation = "rownames")

mg4_c <- runCellQC(mg4_c, sample = NULL,
                 algorithms = c("scDblFinder", "decontX"),
                 geneSetListLocation = "rownames")

wnv_mg3 <- runCellQC(wnv_mg3, sample = NULL,
                 algorithms = c("scDblFinder", "decontX"),
                 geneSetListLocation = "rownames")

wnv_mg4 <- runCellQC(wnv_mg4, sample = NULL,
                 algorithms = c("scDblFinder", "decontX"),
                 geneSetListLocation = "rownames")

#subset:
mg3_c <- subsetSCECols(mg3_c, colData = c("decontX_contamination < 0.6",  
                                      "scDblFinder_doublet_score < 0.9"))

mg4_c <- subsetSCECols(mg4_c, colData = c("decontX_contamination < 0.6",  
                                          "scDblFinder_doublet_score < 0.9"))

wnv_mg3 <- subsetSCECols(wnv_mg3, colData = c("decontX_contamination < 0.6",  
                                          "scDblFinder_doublet_score < 0.9"))

wnv_mg4 <- subsetSCECols(wnv_mg4, colData = c("decontX_contamination < 0.6",  
                                          "scDblFinder_doublet_score < 0.9"))

#convert to Seurat objects:
mg3c <- convertSCEToSeurat(mg3_c)
mg4c <- convertSCEToSeurat(mg4_c)
wnvMg3 <- convertSCEToSeurat(wnv_mg3)
wnvMg4 <- convertSCEToSeurat(wnv_mg4)

#merge, retaining sample ID (add sample condition for later use):
merged.srt <- merge(mg3c, y = c(mg4c, wnvMg3, wnvMg4),
                    add.cell.ids = c("mg3c_mock", "mg4c_mock", "wnvMg3_wnv", "wnvMg4_wnv"),
                    project = "mg_12dpi")

#create 'sample' and 'condition' columns in Seurat object metadata:
merged.srt$sample <- rownames(merged.srt@meta.data)

merged.srt@meta.data <- separate(merged.srt@meta.data, col = 'sample', into = c("sample", "condition", "barcode"),
         sep = "_")

#Seurat QC:
merged.srt[["percent_mt"]] <- PercentageFeatureSet(merged.srt, pattern = "^MT-")
merged.srt <- subset(merged.srt, subset = nFeature_RNA > 100 & 
                       nFeature_RNA < 2500 & 
                       percent_mt < 30)

#explore whether we need to do integration:
merged.srt <- NormalizeData(merged.srt)
merged.srt <- FindVariableFeatures(merged.srt)
merged.srt <- ScaleData(merged.srt)
merged.srt <- RunPCA(merged.srt)
ElbowPlot(merged.srt)
merged.srt <- FindNeighbors(merged.srt, dims = 1:17)
merged.srt <- FindClusters(merged.srt)
merged.srt <- RunUMAP(merged.srt, dims = 1:17)

#save merged sample file for later use:
#saveRDS(merged.srt, file = "/Users/emilyfitzmeyer/Desktop/mito50_merged_mg3mg4.rds", compress = TRUE)

#visualize pre-integrated samples:
sample_plot <- DimPlot(merged.srt, reduction = 'umap', group.by = 'sample')
condition_plot <- DimPlot(merged.srt, reduction = 'umap', group.by = 'condition',
        cols = c('salmon', 'green'))

grid.arrange(sample_plot, condition_plot, ncol = 2)
#(do they appear to be clustering by sample? 
#If so, integrate. If not, proceed to pseudo-bulk w/ merged sample file)

#perform integration:
srt.obj.list <- SplitObject(merged.srt, split.by = 'sample')
for(i in 1:length(srt.obj.list)){
  srt.obj.list[[i]] <- NormalizeData(srt.obj.list[[i]])
  srt.obj.list[[i]] <- FindVariableFeatures(srt.obj.list[[i]])
}

#select integration features:
features <- SelectIntegrationFeatures(srt.obj.list)

#find integration anchors (CCA):
anchors <- FindIntegrationAnchors(srt.obj.list,
                       anchor.features = features)

#integrate data:
srt.integrated <- IntegrateData(anchorset = anchors)

#scale, run PCA, and visualize integrated data:
srt.integrated <- ScaleData(srt.integrated)
srt.integrated <- RunPCA(srt.integrated)
#dims 1:20 is a random selection - I followed a tutorial and she randomly selected 1:50.
srt.integrated <- RunUMAP(srt.integrated, dims = 1:20)

sample_plot2 <- DimPlot(srt.integrated, reduction = 'umap', group.by = 'sample')
condition_plot2 <- DimPlot(srt.integrated, reduction = 'umap', group.by = 'condition',
                          cols = c('salmon', 'green'))

grid.arrange(sample_plot2, condition_plot2, ncol = 2)
grid.arrange(sample_plot, condition_plot, sample_plot2, condition_plot2, ncol = 2, nrow = 2)

#save integrated sample file for later use:
saveRDS(srt.integrated, file = "/Users/emilyfitzmeyer/Desktop/integrated_mg3mg4.rds", compress = TRUE)


#DESeq2 pseudo-bulk:

#srt.obj <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/integrated_mg3mg4.rds")
srt.obj <- readRDS("/Users/emilyfitzmeyer/Desktop/integrated_mg3mg4.rds")

#View(srt.obj@meta.data)

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
saveRDS(dds.results, file = "/Users/emilyfitzmeyer/Desktop/mg3mg4_DESeq2_results.rds", compress = TRUE)
write.csv(dds.results, file = "/Users/emilyfitzmeyer/Desktop/mg3mg4_DESeq2_results.csv")

#dds.results <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq/DE_analysis/pseudo_bulk/mg3mg4_DESeq2_results.rds")

vol1 <- EnhancedVolcano(dds.results, 
                lab = rownames(dds.results),
                x = 'log2FoldChange',
                y = 'pvalue')

