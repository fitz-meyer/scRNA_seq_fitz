library(singleCellTK)
library(tidyverse)
library(Seurat)
library(scCustomize)
library(gridExtra)
library(EnhancedVolcano)
library(clustree)

mg3.c.sce <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/mg3_c_QC.rds")

mg4.c.sce <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/mg4_c_QC.rds")

wnv.mg3.sce <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/wnv_mg3_QC.rds")

wnv.mg4.sce <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/wnv_mg4_QC.rds")

#subset:
mg3.c.sce <- subsetSCECols(mg3.c.sce, colData = c("decontX_contamination < 0.6",  
                                                  "scDblFinder_doublet_score < 0.9"))

mg4.c.sce <- subsetSCECols(mg4.c.sce, colData = c("decontX_contamination < 0.6",  
                                                  "scDblFinder_doublet_score < 0.9"))

wnv.mg3.sce <- subsetSCECols(wnv.mg3.sce, colData = c("decontX_contamination < 0.6",  
                                                      "scDblFinder_doublet_score < 0.9"))

wnv.mg4.sce <- subsetSCECols(wnv.mg4.sce, colData = c("decontX_contamination < 0.6",  
                                                      "scDblFinder_doublet_score < 0.9"))

#convert to Seurat objects and save srt.rds files:
mg3c <- convertSCEToSeurat(mg3.c.sce)
saveRDS(mg3c, "/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/srt_rds/mg3c.rds")
mg4c <- convertSCEToSeurat(mg4.c.sce)
saveRDS(mg4c, "/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/srt_rds/mg4c.rds")
wnvMg3 <- convertSCEToSeurat(wnv.mg3.sce)
saveRDS(wnvMg3, "/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/srt_rds/wnvMg3.rds")
wnvMg4 <- convertSCEToSeurat(wnv.mg4.sce)
saveRDS(wnvMg4, "/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/srt_rds/wnvMg4.rds")

#read srt RDS files in if restarting srt workflow
mg3c <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/srt_rds/mg3c.rds")
mg4c <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/srt_rds/mg4c.rds")
wnvMg3 <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/srt_rds/wnvMg3.rds")
wnvMg4 <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/srt_rds/wnvMg4.rds")

#merge, retaining sample ID (add sample condition for later use):
merged.srt <- merge(mg3c, y = c(mg4c, wnvMg3, wnvMg4),
                    add.cell.ids = c("mg3c_mock", "mg4c_mock", "wnvMg3_wnv", "wnvMg4_wnv"),
                    project = "mg_12dpi")

merged.srt$sample <- rownames(merged.srt@meta.data)

merged.srt@meta.data <- separate(merged.srt@meta.data, col = 'sample', into = c("sample", "condition", "barcode"),
                                 sep = "_")

#Seurat QC:
merged.srt[["percent_mt"]] <- PercentageFeatureSet(merged.srt, pattern = "^MT-")
merged.srt <- subset(merged.srt, subset = nFeature_RNA > 100 & 
                       nFeature_RNA < 2500 & 
                       percent_mt < 25)

#generate and transpose counts matrix for scLink:
# split.srt1 <- SplitObject(merged.srt, split.by = "condition")
# wnv.srt1 <- split.srt1$wnv
# unique(wnv.srt1$sample)
# counts.matrix <- wnv.srt1@assays$RNA@counts
# t.12dpi.counts.matrix <- t(counts.matrix)
# saveRDS(t.12dpi.counts.matrix, "/Users/emilyfitzmeyer/Desktop/12dpi_wnv_counts_matrix.rds")

merged.srt <- NormalizeData(merged.srt)
merged.srt <- FindVariableFeatures(merged.srt)

# wnv.srt1 <- NormalizeData(wnv.srt1)
# wnv.srt1 <- FindVariableFeatures(wnv.srt1)

#making gene of interest list for scLink:
# vf.wnv12dpi.srt <- VariableFeatures(wnv.srt1)
# vf.12dpi.top500 <- head(vf.wnv12dpi.srt, 500)
# saveRDS(vf.12dpi.top500, "/Users/emilyfitzmeyer/Desktop/vf_wnv12dpi_top500.rds")

#scale:
#srt_genes <- rownames(srt)
#optional "features = srt_genes" argument
merged.srt <- ScaleData(merged.srt)
merged.srt <- RunPCA(merged.srt)
#ElbowPlot(merged.srt, ndims = 40)
merged.srt <- FindNeighbors(merged.srt, dims = 1:35)
# merged.srt <- FindClusters(merged.srt, resolution = c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.6))
# clustree(merged.srt)
# setwd("/Users/emilyfitzmeyer/Desktop/")
# ggsave("clustree_12dpi.png", plot = last_plot(), device = png(), scale = 1, width = 9.5, height = 9.5, dpi = 300)
# dev.off()
merged.srt <- FindClusters(merged.srt, resolution = 0.6)
merged.srt <- RunUMAP(merged.srt, dims = 1:35)
#save merged object:
saveRDS(merged.srt, file = "/Users/emilyfitzmeyer/Desktop/12dpi_res0.6.rds")
#Read in merged.srt
merged.srt <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/merged_rds/12dpi_res1.1.rds")
#for terminal
#merged.srt <- readRDS("/home/fitzmeyer/rds/12dpi_res1.1.rds")

#rename clusters
new_idents <- c("EC-like-1", "1", "2", "EC-like-2", "4", "ISC/EB", "EC-like-3", "EC", "8", "VM", "10", "cardia", "12", "EE", "ISC/EB-prol", "15", "HC")
merged.srt <- Rename_Clusters(merged.srt, new_idents = new_idents, meta_col_name = "res1.1_idents")

#Examine clusters in each condition
clusters_12dpi <- DimPlot(merged.srt, reduction = "umap", label = TRUE, repel = TRUE, split.by = "condition", pt.size = 0.2) +
  theme(axis.text.x = element_text(size = 15)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(legend.position = 'none')

ggsave("clusters_12dpi.png", plot = last_plot(), device = png(), scale = 1, width = 8.5, height = 5, dpi = 300)
dev.off()

#Examine clusters in each condition - NO SPLIT
clusters_12dpi <- DimPlot(merged.srt, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.2) +
  theme(axis.text.x = element_text(size = 15)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(legend.position = 'none')

ggsave("clusters_12dpi.png", plot = last_plot(), device = png(), scale = 1, width = 4.25, height = 5, dpi = 300)
dev.off()

#CLUSTER MARKER DOTPLOT
dotplot_genes <- c("gene11957", "gene316", "gene769", "gene10632", "gene2985", "gene9804", "gene11104", "gene2439", 
                   "gene452", "gene7512", "gene1416", "gene13104", "gene2167", "gene12446", "gene9049", "gene3360", "gene3283", "gene11820", "gene6403", 
                   "gene6861", "gene13589", "gene14180", "gene5110", "gene4641", "gene12764",
                   "gene10257", "gene927")

dotplot_break_labels <- c("POU2F1", "PLA2G6", "AGBL5", "PROX1", "ACTB", "Mlc2", "Mhc", "NIMB2", "SPARC", "pebIII_CPIJ002629", "pebIII_CPIJ002609", 
                          "klumpfuss", "PCNA", "C-type lysozyme", "Sugar_tr_CPIJ019592", "Sugar_tr_CPIJ014327", "Sugar_tr_CPIJ011910",
                          "Sugar_tr_CPIJ012675", "Sugar_tr_CPIJ012678", "chitin-binding_CPIJ004734", "chitin-binding_CPIJ004734",
                          "serine_protease_CPIJ006568", "serine_protease_CPIJ015103",
                          "serine_protease_CPIJ007079", "serine_protease_KDR18614", 
                          "Mal-B2", "Mal-A4")

DotPlot_scCustom(merged.srt, features = dotplot_genes) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = dotplot_break_labels)

setwd("/Users/emilyfitzmeyer/Desktop/")
ggsave("cl_marker_dotplot_12dpi.png", plot = last_plot(), device = png(), scale = 1, width = 9, height = 5.5, dpi = 300)
dev.off()

#cl3vcl10 dotplot
dotplot_genes <- c("nbisL1-trna-9", "MT-nbis-gene-3", "MT-nbis-gene-2", "nbisL1-trna-7", "nbisL1-trna-10",
                   "nbisL1-trna-6", "nbisL1-trna-16", "nbisL1-trna-3", "nbisL1-trna-17", "MT-nbis-gene-5",
                   "nbisL1-trna-20", "nbisL1-trna-19", "nbis-gene-2-utr")

dotplot_break_labels <- c("MT-tRNA-Asp", "COX3", "ATP6", "MT-tRNA-Leu", "MT-tRNA-Gly",
                          "MT-tRNA-Tyr", "MT-tRNA-Phe", "MT-tRNA-Met", "MT-tRNA-His", "CYTB",
                          "MT-tRNA-Ser", "MT-tRNA-Pro", "WNV 5' UTR")

DotPlot_scCustom(split.srt$wnv, features = dotplot_genes) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = dotplot_break_labels)

ggsave("cl10_marker_dotplot.png", plot = last_plot(), device = png(), scale = 1, width = 6, height = 4, dpi = 300)
dev.off()

#ISC/EB clusters dotplot
dotplot_genes <- c("gene6060", "gene4986", "gene10089",	"gene6015",	"gene6292",	"gene1252",	"gene6564",
                   "gene2793", "gene6552", "gene5534", "gene4419")

dotplot_break_labels <- c("DUF4803", "gene4986", "cecropin", "cecropin", "CUTA", "GSTD11", "awd",
                          "NHP2L1", "histone H2A", "RAN", "BIRC5")

DotPlot_scCustom(merged.srt, features = dotplot_genes) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = dotplot_break_labels)
ggsave("ISC_marker_dotplot.png", plot = last_plot(), device = png(), scale = 1, width = 6, height = 4, dpi = 300)
dev.off()

#show sample and condition induced impact on cellular heterogeneity:
sample_plot <- DimPlot(merged.srt, reduction = 'umap', group.by = 'sample') &
  theme(axis.text.x = element_text(size = 17)) &
  theme(axis.text.y = element_text(size = 17)) &
  ggtitle("12dpi sample grouping")
ggsave("sample_grouping_12dpi.png", plot = last_plot(), device = png(), scale = 1, width = 6.5, height = 5, dpi = 300)
dev.off()

condition_plot <- DimPlot(merged.srt, reduction = 'umap', group.by = 'condition',
                          cols = c('deepskyblue', 'salmon')) &
  theme(axis.text.x = element_text(size = 17)) &
  theme(axis.text.y = element_text(size = 17)) &
  ggtitle("12dpi condition grouping")
ggsave("condition_grouping_12dpi.png", plot = last_plot(), device = png(), scale = 1, width = 6.5, height = 5, dpi = 300)
dev.off()

#split obj for feature plots and marker identification
split.srt <- SplitObject(merged.srt, split.by = "condition")

mean(split.srt$mock$nFeature_RNA)

#plot presence of WNV 5' UTR:
WNV.utr.plot <- FeaturePlot(merged.srt, features = "nbis-gene-2-utr", split.by = "condition", pt.size = 0.5) &
  theme(axis.text.x = element_text(size = 17)) &
  theme(axis.text.y = element_text(size = 17)) &
  theme(legend.position = c(0.9, 0.22))

setwd("/Users/emilyfitzmeyer/Desktop/")
ggsave("wnv_utr_12dpi.png", plot = last_plot(), device = png(), scale = 1, width = 10, height = 5, dpi = 300)
dev.off()

mock.percent.whole <- Percent_Expressing(split.srt$mock, features = "nbis-gene-2-utr", entire_object = TRUE)
wnv.percent.whole <- Percent_Expressing(split.srt$wnv, features = "nbis-gene-2-utr", entire_object = TRUE)
wnv.percent.cluster <- Percent_Expressing(split.srt$wnv, features = "nbis-gene-2-utr")
wnv.expLevel.whole <- AverageExpression(split.srt$wnv, features = "nbis-gene-2-utr", group.by = "orig.ident")


VlnPlot(split.srt$wnv, features = "nbis-gene-2-utr", cols = c("darkslategray3","darkslategray3","darkslategray3","darkslategray3","darkslategray3","darkslategray3","darkslategray3","darkslategray3","darkslategray3","darkslategray3","darkslategray3","darkslategray3","darkslategray3","darkslategray3","darkslategray3","darkslategray3","darkslategray3")) + 
  theme(axis.text.x = element_text(size = 15)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(legend.position = 'none') +
  scale_x_discrete(labels = new_idents) +
  labs(x = NULL)
ggsave("12dpi_wnv_pct_cluster.png", plot = last_plot(), device = png(), scale = 1, width = 5.5, height = 3.5, dpi = 300)
dev.off()

plot_temp <- VlnPlot(split.srt$wnv, features = "nbis-gene-2-utr", idents = c("0", "3", "5", "14", "6", "7", "9", "11", "13", "16")) + 
  # ylim(0, 20) +
  # stat_compare_means(comparisons = my_comparisons, label = "p.signif", hide.ns = TRUE) +
  theme(legend.position = 'none')

write.csv(plot_temp$data, file = "/Users/emilyfitzmeyer/Desktop/knwn_CLs_WNV_vRNA_12dpi.csv")

ggsave("12dpi_wnv_pct_knwnCluster.png", plot = last_plot(), device = png(), scale = 1, width = 6, height = 6, dpi = 300)
dev.off()


#subset clusters:
ISC.srt <- subset(merged.srt, idents = c("ISC/EB", "ISC/EB-prol"))

#VISUALIZE SPECIFIC GENES
FeatureScatter(split.srt$wnv, feature1 = "nbis-gene-2-utr", feature2 = "gene11340", plot.cor = FALSE) + 
  theme(legend.position = "none") +
  labs(y = "IMD", x = "WNV 5' UTR")

VlnPlot(merged.srt, features = "gene11340", group.by = "condition") +
  theme(legend.position = "none") +
  ggtitle("IMD")

#ISC plots
FeaturePlot(merged.srt, features = "gene12926", pt.size = 0.5) +
  theme(axis.text.x = element_text(size = 17)) +
  theme(axis.text.y = element_text(size = 17)) +
  scale_x_continuous(limits = c(-12, -6)) +
  scale_y_continuous(limits = c(-2, 3)) +
  theme(legend.position = c(0.9, 0.23)) +
  theme(legend.text = element_text(size = 17)) +
  ggtitle("AURKB")
ggsave("AURKB_12dpi.png", plot = last_plot(), device = png(), scale = 1, width = 6, height = 5, dpi = 300)
dev.off()

#PROX1 plots
FeaturePlot(merged.srt, features = "gene10632") +
  theme(axis.text.x = element_text(size = 17)) +
  theme(axis.text.y = element_text(size = 17)) +
  ggtitle("PROX1")
ggsave("PROX1_12dpi.png", plot = last_plot(), device = png(), scale = 1, width = 6, height = 6, dpi = 300)
dev.off()
VlnPlot(merged.srt, features = "gene10632") + 
  theme(legend.position = 'none') + 
  theme(axis.text.x = element_text(size = 17)) +
  theme(axis.text.y = element_text(size = 17)) +
  ggtitle("PROX1")
ggsave("PROX1_vlnPlot_12dpi.png", plot = last_plot(), device = png(), scale = 1, width = 6, height = 6, dpi = 300)
dev.off()

#CLUSTER PROPORTIONS
cluster_proportions <- prop.table(table(Idents(merged.srt)))
as.data.frame(cluster_proportions)




#identifying pctExp>75 genes for COG analysis
#my adaptation - chatGPT version wasn't quite doing what I wanted
#note - my version needs to be run on the server bc it's chonky 
clusters <- unique(merged.srt$seurat_clusters)
gene_matrices <- list()

for (cluster in clusters) {
  subset_seurat <- subset(merged.srt, ident = cluster)
  cluster_genes <- rownames(subset_seurat@assays$RNA@data)
  pct_exp <- Percent_Expressing(subset_seurat, features = cluster_genes, entire_object = TRUE)
  top_80pct_genes <- subset(pct_exp, pct_exp$All_Cells >= 75)
  gene_matrices[[cluster]] <- top_80pct_genes
  print("next cluster")
}

saveRDS(gene_matrices, file = "/home/fitzmeyer/rds/gene_matrices_12dpi.rds")

cl_nums <- names(gene_matrices)
cl_nums <- paste0(cl_nums, "pctExp80.csv")
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


#percent expressing key immune genes:
#===============================================================================
immn.gene.pctExp <- Percent_Expressing(merged.srt, assay = 'RNA', split_by = "condition", 
                                       features = c("gene1501", "gene7624", "gene10322", "gene11340", "gene3044",
                                                    "gene6258", "gene8444", "gene11364", "gene6091", "gene6092",
                                                    "gene5869", "gene5870", "gene10292", "gene10293", "gene3397",
                                                    "gene14111", "gene2893", "gene3913", "gene2450", "gene5218", 
                                                    "gene13590", "gene3988", "gene3989", "gene4535", "gene14091"),
                                       entire_object = TRUE)
  
#averageExpression key immune genes:
#===============================================================================

immn.gene.avgExp <- AverageExpression(merged.srt, assays = 'RNA', 
                                      features = c("gene1501", "gene7624", "gene10322", "gene11340", "gene3044",
                                                   "gene6258", "gene8444", "gene11364", "gene6091", "gene6092",
                                                   "gene5869", "gene5870", "gene10292", "gene10293", "gene3397",
                                                   "gene14111", "gene2893", "gene3913", "gene2450", "gene5218", 
                                                   "gene13590", "gene3988", "gene3989", "gene4535", "gene14091"), 
                                      group.by = "condition")

#CLUSTER MARKERS
#===============================================================================
#set default assay to "RNA" so DE analyses use raw counts:
DefaultAssay(merged.srt) <- "RNA"
split.srt <- SplitObject(merged.srt, split.by = "condition")
#Find all cluster markers for wnv+ samples:
markers_all <- FindAllMarkers(merged.srt, 
                          logfc.threshold = 0.25,
                          test.use = "wilcox")

write.csv(markers_all, file = "/Users/emilyfitzmeyer/Desktop/12dpi_all_markers_wilcox.csv")

#and mock samples:
markers_mock <- FindAllMarkers(split.srt$mock, 
                          logfc.threshold = 0.25,
                          test.use = "wilcox")

write.csv(markers_mock, file = "/Users/emilyfitzmeyer/Desktop/mock_12dpi_all_markers_wilcox.csv")

#identify all/conserved cell type markers:
cnsrv.markers <- FindConservedMarkers(merged.srt, 
                                      ident.1 = 16,
                                      grouping.var = "condition",
                                      slot = "data",
                                      test.use = "wilcox",
                                      verbose = FALSE)

write.csv(cnsrv.markers, file = "/Users/emilyfitzmeyer/Desktop/cl16_conservedMarkers_wilcox.csv")

#find markers that differentiate clusters 
cl_markers <- FindMarkers(split.srt$wnv,
                                 ident.1 = 1,
                                 slot = "data",
                                 test.use = "wilcox",
                                 verbose = FALSE)

write.csv(cl_markers, file = "/Users/emilyfitzmeyer/Desktop/cl1_v_all_markers_WNVonly_wilcox.csv")



#DESeq2 pseudo-bulk 
#===============================================================================
#srt.obj <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/integrated_mg3mg4.rds")
srt.obj <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/merged_rds/12dpi_res0.6.rds")

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
saveRDS(dds.results, file = "/Users/emilyfitzmeyer/Desktop/12dpi_DESeq2_results.rds", compress = TRUE)
write.csv(dds.results, file = "/Users/emilyfitzmeyer/Desktop/12dpi_DESeq2_results.csv")

dds.results <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/merged_rds/12dpi_DESeq2_results.rds")
results <- read.csv("/Users/emilyfitzmeyer/Desktop/scRNAseq/plots/12dpi_labChange_volcanoPlotValues.csv")

#cutoff varies based on results file
EnhancedVolcano(results, 
                lab = results$gene_ID,
                title = "WNV vs. Mock DEGs 12dpi",
                subtitle = NULL,
                caption = NULL,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.000065,
                drawConnectors = TRUE,
                arrowheads = FALSE,
                labSize = 5,
                axisLabSize = 23,
                legendIconSize = 3)

setwd("/Users/emilyfitzmeyer/Desktop/")
ggsave("PB_DEGs_12dpi.png", plot = last_plot(), device = png(), scale = 1, width = 6, height = 6, dpi = 300)
dev.off()



