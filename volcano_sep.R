library(singleCellTK)
library(tidyverse)
library(Seurat)
library(scCustomize)
library(gridExtra)
library(EnhancedVolcano)
library(clustree)

results <- read.csv("/Users/emilyfitzmeyer/Desktop/scRNAseq/figures/WNV_specific_12v4_DEG_values.csv")

#cutoff varies based on results file
EnhancedVolcano(results, 
                lab = results$gene_ID,
                title = "12dpi vs. 4dpi infection specific DEGs",
                subtitle = NULL,
                caption = NULL,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.001,
                FCcutoff = 2.0,
                drawConnectors = TRUE,
                arrowheads = FALSE,
                labSize = 5,
                axisLabSize = 23,
                legendIconSize = 3)

setwd("/Users/emilyfitzmeyer/Desktop/")
ggsave("PB_DEGs_12dpi.png", plot = last_plot(), device = png(), scale = 1, width = 6, height = 6, dpi = 300)
dev.off()




