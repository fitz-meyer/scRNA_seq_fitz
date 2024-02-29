library(tidyverse)
library(Seurat)
library(gridExtra)
library(ggplot2)


merged.srt <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/merged_rds/4dpi_res0.6.rds")
HC.srt <- subset(merged.srt, idents = 12)

#make binary version of "data" slot
bin_dat_HC <- GetAssayData(HC.srt[["RNA"]], slot = "data")
bin_dat_HC[bin_dat_HC>0] <- 1
HC.srt[["binary"]] <- CreateAssayObject(data = bin_dat_HC)
DefaultAssay(HC.srt) <- "binary"

#Mature granulocyte:

gene9384_4dpi_12 <- FeaturePlot(HC.srt, features = "gene9384", pt.size = 0.5, cols = c("grey", "grey")) +
  theme(axis.text.x = element_text(size = 15, angle = 45, vjust = 0.5)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(legend.position = 'none') +
  ggtitle("SCRASP1") +
  scale_x_continuous(limits = c(7.5, 9)) +
  scale_y_continuous(limits = c(0, 2)) +
  labs(x = NULL, y = "HC-1")

#Oenocytoid

#gene6205 <- FeaturePlot(HC.srt, features = "gene6205", pt.size = 0.5, cols = c("lightgrey", "lightgrey")) &
gene6205_4dpi_12 <- FeaturePlot(HC.srt, features = "gene6205", pt.size = 0.5) +
  theme(axis.text.x = element_text(size = 15, angle = 45, vjust = 0.5)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(legend.position = 'none') +
  ggtitle("SCRB3") +
  scale_x_continuous(limits = c(7.5, 9)) +
  scale_y_continuous(limits = c(0, 2), labels = NULL) +
  labs(x = NULL, y = NULL)

#Universal:

gene2439_4dpi_12 <- FeaturePlot(HC.srt, features = "gene2439", pt.size = 0.5, cols = c("grey", "grey")) +
  theme(axis.text.x = element_text(size = 15, angle = 45, vjust = 0.5)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(legend.position = 'none') +
  ggtitle("NIMB2") +
  scale_x_continuous(limits = c(7.5, 9)) +
  scale_y_continuous(limits = c(0, 2), labels = NULL) +
  labs(x = NULL, y = NULL)

gene452_4dpi_12 <- FeaturePlot(HC.srt, features = "gene452", pt.size = 0.5) +
  theme(axis.text.x = element_text(size = 15, angle = 45, vjust = 0.5)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(legend.position = 'none') +
  ggtitle("SPARC") +
  scale_x_continuous(limits = c(7.5, 9)) +
  scale_y_continuous(limits = c(0, 2), labels = NULL) +
  labs(x = NULL, y = NULL)

#NEW CLUSTER

HC_2.srt <- subset(merged.srt, idents = 13)

#make binary version of "data" slot
bin_dat_HC_2 <- GetAssayData(HC_2.srt[["RNA"]], slot = "data")
bin_dat_HC_2[bin_dat_HC_2>0] <- 1
HC_2.srt[["binary"]] <- CreateAssayObject(data = bin_dat_HC_2)
DefaultAssay(HC_2.srt) <- "binary"

#Mature granulocyte:

gene9384_4dpi_13 <- FeaturePlot(HC_2.srt, features = "gene9384", pt.size = 0.5) +
  theme(axis.text.x = element_text(size = 15, angle = 45, vjust = 0.5)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(legend.position = 'none') +
  ggtitle(NULL) +
  scale_x_continuous(limits = c(13, 14.5)) +
  scale_y_continuous(limits = c(2, 4)) +
  labs(x = NULL, y = "HC-2")

#Oenocytoid

gene6205_4dpi_13 <- FeaturePlot(HC_2.srt, features = "gene6205", pt.size = 0.5, cols = c("grey", "grey")) +
  theme(axis.text.x = element_text(size = 15, angle = 45, vjust = 0.5)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(legend.position = 'none') +
  ggtitle(NULL) +
  scale_x_continuous(limits = c(13, 14.5)) +
  scale_y_continuous(limits = c(2, 4), labels = NULL) +
  labs(x = NULL, y = NULL)

#Universal:
gene2439_4dpi_13 <- FeaturePlot(HC_2.srt, features = "gene2439", pt.size = 0.5, cols = c("blue", "blue")) +
  theme(axis.text.x = element_text(size = 15, angle = 45, vjust = 0.5)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(legend.position = 'none') +
  ggtitle(NULL) +
  scale_x_continuous(limits = c(13, 14.5)) +
  scale_y_continuous(limits = c(2, 4), labels = NULL) +
  labs(x = NULL, y = NULL)

gene452_4dpi_13 <- FeaturePlot(HC_2.srt, features = "gene452", pt.size = 0.5, cols = c("blue", "blue")) +
  theme(axis.text.x = element_text(size = 15, angle = 45, vjust = 0.5)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(legend.position = 'none') +
  ggtitle(NULL) +
  scale_x_continuous(limits = c(13, 14.5)) +
  scale_y_continuous(limits = c(2, 4), labels = NULL) +
  labs(x = NULL, y = NULL)

#NEW CLUSTER
merged.srt2 <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq/rds_files/merged_rds/12dpi_res1.1.rds")
HC_3.srt <- subset(merged.srt2, idents = 16)

#make binary version of "data" slot
bin_dat_HC_3 <- GetAssayData(HC_3.srt[["RNA"]], slot = "data")
bin_dat_HC_3[bin_dat_HC_3>0] <- 1
HC_3.srt[["binary"]] <- CreateAssayObject(data = bin_dat_HC_3)
DefaultAssay(HC_3.srt) <- "binary"

#Mature granulocyte:

gene9384_12dpi_16 <- FeaturePlot(HC_3.srt, features = "gene9384", pt.size = 0.5) +
  theme(axis.text.x = element_text(size = 15, angle = 45, vjust = 0.5)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(legend.position = 'none') +
  ggtitle(NULL) +
  scale_x_continuous(limits = c(-2, 1)) +
  scale_y_continuous(limits = c(2, 8)) +
  labs(x = NULL, y = "HC")

#Oenocytoid

gene6205_12dpi_16 <- FeaturePlot(HC_3.srt, features = "gene6205", pt.size = 0.5) +
  theme(axis.text.x = element_text(size = 15, angle = 45, vjust = 0.5)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(legend.position = 'none') +
  ggtitle(NULL) +
  scale_x_continuous(limits = c(-2, 1)) +
  scale_y_continuous(limits = c(2, 8), labels = NULL) +
  labs(x = NULL, y = NULL)

#Universal:

gene2439_12dpi_16 <- FeaturePlot(HC_3.srt, features = "gene2439", pt.size = 0.5) +
  theme(axis.text.x = element_text(size = 15, angle = 45, vjust = 0.5)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(legend.position = 'none') +
  ggtitle(NULL) +
  scale_x_continuous(limits = c(-2, 1)) +
  scale_y_continuous(limits = c(2, 8), labels = NULL) +
  labs(x = NULL, y = NULL)

gene452_12dpi_16 <- FeaturePlot(HC_3.srt, features = "gene452", pt.size = 0.5) +
  theme(axis.text.x = element_text(size = 15, angle = 45, vjust = 0.5)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(legend.position = 'none') +
  ggtitle(NULL) +
  scale_x_continuous(limits = c(-2, 1)) +
  scale_y_continuous(limits = c(2, 8), labels = NULL) +
  labs(x = NULL, y = NULL)


#COMBINE ALL PLOTS

plot <- gene9384_4dpi_12 + gene6205_4dpi_12 + gene2439_4dpi_12 + gene452_4dpi_12 + gene9384_4dpi_13 + gene6205_4dpi_13 + gene2439_4dpi_13 + gene452_4dpi_13 + gene9384_12dpi_16 + gene6205_12dpi_16 + gene2439_12dpi_16 + gene452_12dpi_16 + patchwork::plot_layout(ncol = 4)

setwd("/Users/emilyfitzmeyer/Desktop/")
ggsave("HC_plot.png", plot = plot, device = png(), scale = 1, width = 9.4, height = 7, dpi = 300)
dev.off()




plot2 <- clusters_4dpi + clusters_12dpi + patchwork::plot_layout(ncol = 1)
ggsave("clusters.png", plot = plot2, device = png(), scale = 1, width = 8.5, height = 10, dpi = 300)

