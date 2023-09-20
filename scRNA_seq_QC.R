library(singleCellTK)
library(tidyverse)
library(Seurat)

r10x <- Read10X(data.dir = "/Users/emilyfitzmeyer/Desktop/scRNAseq/cellranger_aggr/20230916_AGGR_mg3mg4_count/outs/filtered_feature_bc_matrix/")
srt.int <- CreateSeuratObject(counts = r10x)

sce <- importCellRanger(sampleDirs = "/Users/emilyfitzmeyer/Desktop/scRNAseq/cellranger_aggr/20230916_AGGR_mg3mg4_count/",
                         dataType = c("filtered"),
                         matrixFileNames = "matrix.mtx.gz",
                         featuresFileNames = "features.tsv.gz",
                         barcodesFileNames = "barcodes.tsv.gz",
                         gzipped = TRUE)

#mito_genes <- c("MT-nbis-gene-5", "MT-nbisL1-trna-20", "MT-nbisL1-trna-6", 
#                "MT-nbisL1-trna-3", "MT-nbisL1-trna-7", "MT-nbisL1-trna-9", 
#                "MT-nbis-gene-2", "MT-nbis-gene-3", "MT-nbisL1-trna-10", 
#                "MT-nbisL1-trna-16", "MT-nbisL1-trna-17", "MT-nbis-gene-4", 
#                "MT-nbisL1-trna-19")

#sce <- runCellQC(sce, algorithms = "QCMetrics", mitoGeneLocation = NULL, mitoPrefix="$$", geneSetList = list(mito = mito_genes), geneSetListLocation = "rownames")

sce <- runCellQC(sce, sample = NULL,
                 algorithms = c("scDblFinder", "decontX"),
                 geneSetListLocation = "rownames")

saveRDS(sce, file = "Desktop/scRNAseq/QC_rds_files/aggr_mg3mg4_QC.rds", compress = TRUE)

