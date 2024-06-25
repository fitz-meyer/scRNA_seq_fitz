library(tidyverse)
library(Seurat)
library(gridExtra)
library(ggplot2)


merged.srt <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq_pub/mergeDPI_oldNorm/alldpi_oldNorm_nmn_res0.6.rds")
split.srt <- SplitObject(merged.srt, split.by = "dpi")

#NEW DPI

dpi.srt <- split.srt$"4dpi"
HC.srt <- subset(dpi.srt, idents = 14)

#make binary version of "data" slot
bin_dat_HC <- GetAssayData(HC.srt[["RNA"]], slot = "data")
bin_dat_HC[bin_dat_HC>0] <- 1
HC.srt[["binary"]] <- CreateAssayObject(data = bin_dat_HC)
DefaultAssay(HC.srt) <- "binary"

#Mature granulocyte:

gene9384_4dpi_14 <- FeaturePlot(HC.srt, features = "gene9384", pt.size = 0.5) +
  theme(axis.text.x = element_text(size = 13, angle = 45, vjust = 0.5)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(legend.position = 'none') +
  ggtitle("SCRASP1") +
  scale_x_continuous(limits = c(-13.25, -12.3)) +
  scale_y_continuous(limits = c(-0.65, 0.1)) +
  labs(x = NULL, y = "HC-1")

#Oenocytoid

gene6205_4dpi_14 <- FeaturePlot(HC.srt, features = "gene6205", pt.size = 0.5, cols = c("grey", "grey")) +
  theme(axis.text.x = element_text(size = 13, angle = 45, vjust = 0.5)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(legend.position = 'none') +
  ggtitle("SCRB3") +
  scale_x_continuous(limits = c(-13.25, -12.3)) +
  scale_y_continuous(limits = c(-0.65, 0.1), labels = NULL) +
  labs(x = NULL, y = NULL)

#Universal:

gene2439_4dpi_14 <- FeaturePlot(HC.srt, features = "gene2439", pt.size = 0.5, cols = c("blue", "blue")) +
  theme(axis.text.x = element_text(size = 13, angle = 45, vjust = 0.5)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(legend.position = 'none') +
  ggtitle("NIMB2") +
  scale_x_continuous(limits = c(-13.25, -12.3)) +
  scale_y_continuous(limits = c(-0.65, 0.1), labels = NULL) +
  labs(x = NULL, y = NULL)

gene452_4dpi_14 <- FeaturePlot(HC.srt, features = "gene452", pt.size = 0.5, cols = c("blue", "blue")) +
  theme(axis.text.x = element_text(size = 13, angle = 45, vjust = 0.5)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(legend.position = 'none') +
  ggtitle("SPARC") +
  scale_x_continuous(limits = c(-13.25, -12.3)) +
  scale_y_continuous(limits = c(-0.65, 0.1), labels = NULL) +
  labs(x = NULL, y = NULL)


#NEW CLUSTER
HC_2.srt <- subset(dpi.srt, idents = 16)

#make binary version of "data" slot
bin_dat_HC_2 <- GetAssayData(HC_2.srt[["RNA"]], slot = "data")
bin_dat_HC_2[bin_dat_HC_2>0] <- 1
HC_2.srt[["binary"]] <- CreateAssayObject(data = bin_dat_HC_2)
DefaultAssay(HC_2.srt) <- "binary"

#Mature granulocyte:

gene9384_4dpi_16 <- FeaturePlot(HC_2.srt, features = "gene9384", pt.size = 0.5) +
  theme(axis.text.x = element_text(size = 13, angle = 45, vjust = 0.5)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(legend.position = 'none') +
  ggtitle(NULL) +
  scale_x_continuous(limits = c(-7.5, -6)) +
  scale_y_continuous(limits = c(4.5, 6.5)) +
  labs(x = NULL, y = "HC-2")

#Oenocytoid

gene6205_4dpi_16 <- FeaturePlot(HC_2.srt, features = "gene6205", pt.size = 0.5) +
  theme(axis.text.x = element_text(size = 13, angle = 45, vjust = 0.5)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(legend.position = 'none') +
  ggtitle(NULL) +
  scale_x_continuous(limits = c(-7.5, -6)) +
  scale_y_continuous(limits = c(4.5, 6.5), labels = NULL) +
  labs(x = NULL, y = NULL)

#Universal:
gene2439_4dpi_16 <- FeaturePlot(HC_2.srt, features = "gene2439", pt.size = 0.5, cols = c("grey", "grey")) +
  theme(axis.text.x = element_text(size = 13, angle = 45, vjust = 0.5)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(legend.position = 'none') +
  ggtitle(NULL) +
  scale_x_continuous(limits = c(-7.5, -6)) +
  scale_y_continuous(limits = c(4.5, 6.5), labels = NULL) +
  labs(x = NULL, y = NULL)

gene452_4dpi_16 <- FeaturePlot(HC_2.srt, features = "gene452", pt.size = 0.5) +
  theme(axis.text.x = element_text(size = 13, angle = 45, vjust = 0.5)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(legend.position = 'none') +
  ggtitle(NULL) +
  scale_x_continuous(limits = c(-7.5, -6)) +
  scale_y_continuous(limits = c(4.5, 6.5), labels = NULL) +
  labs(x = NULL, y = NULL)

#NEW DPI

dpi.srt <- split.srt$"12dpi"
HC.srt <- subset(dpi.srt, idents = 14)

#make binary version of "data" slot
bin_dat_HC <- GetAssayData(HC.srt[["RNA"]], slot = "data")
bin_dat_HC[bin_dat_HC>0] <- 1
HC.srt[["binary"]] <- CreateAssayObject(data = bin_dat_HC)
DefaultAssay(HC.srt) <- "binary"

#Mature granulocyte:

gene9384_12dpi_14 <- FeaturePlot(HC.srt, features = "gene9384", pt.size = 0.5) +
  theme(axis.text.x = element_text(size = 13, angle = 45, vjust = 0.5)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(legend.position = 'none') +
  ggtitle("SCRASP1") +
  scale_x_continuous(limits = c(-13.25, -12.3)) +
  scale_y_continuous(limits = c(-0.65, 0.1)) +
  labs(x = NULL, y = "HC-1")

#Oenocytoid

gene6205_12dpi_14 <- FeaturePlot(HC.srt, features = "gene6205", pt.size = 0.5, cols = c("grey", "grey")) +
#gene6205_12dpi_14 <- FeaturePlot(HC.srt, features = "gene6205", pt.size = 0.5) +
  theme(axis.text.x = element_text(size = 13, angle = 45, vjust = 0.5)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(legend.position = 'none') +
  ggtitle("SCRB3") +
  scale_x_continuous(limits = c(-13.25, -12.3)) +
  scale_y_continuous(limits = c(-0.65, 0.1), labels = NULL) +
  labs(x = NULL, y = NULL)

#Universal:

gene2439_12dpi_14 <- FeaturePlot(HC.srt, features = "gene2439", pt.size = 0.5) +
  theme(axis.text.x = element_text(size = 13, angle = 45, vjust = 0.5)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(legend.position = 'none') +
  ggtitle("NIMB2") +
  scale_x_continuous(limits = c(-13.25, -12.3)) +
  scale_y_continuous(limits = c(-0.65, 0.1), labels = NULL) +
  labs(x = NULL, y = NULL)

gene452_12dpi_14 <- FeaturePlot(HC.srt, features = "gene452", pt.size = 0.5, cols = c("blue", "blue")) +
  theme(axis.text.x = element_text(size = 13, angle = 45, vjust = 0.5)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(legend.position = 'none') +
  ggtitle("SPARC") +
  scale_x_continuous(limits = c(-13.25, -12.3)) +
  scale_y_continuous(limits = c(-0.65, 0.1), labels = NULL) +
  labs(x = NULL, y = NULL)


#NEW CLUSTER
HC_2.srt <- subset(dpi.srt, idents = 16)

#make binary version of "data" slot
bin_dat_HC_2 <- GetAssayData(HC_2.srt[["RNA"]], slot = "data")
bin_dat_HC_2[bin_dat_HC_2>0] <- 1
HC_2.srt[["binary"]] <- CreateAssayObject(data = bin_dat_HC_2)
DefaultAssay(HC_2.srt) <- "binary"

#Mature granulocyte:

gene9384_12dpi_16 <- FeaturePlot(HC_2.srt, features = "gene9384", pt.size = 0.5, cols = c("grey", "grey")) +
  theme(axis.text.x = element_text(size = 13, angle = 45, vjust = 0.5)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(legend.position = 'none') +
  ggtitle(NULL) +
  scale_x_continuous(limits = c(-7.5, -6)) +
  scale_y_continuous(limits = c(4.5, 6.5)) +
  labs(x = NULL, y = "HC-2")

#Oenocytoid

gene6205_12dpi_16 <- FeaturePlot(HC_2.srt, features = "gene6205", pt.size = 0.5) +
  theme(axis.text.x = element_text(size = 13, angle = 45, vjust = 0.5)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(legend.position = 'none') +
  ggtitle(NULL) +
  scale_x_continuous(limits = c(-7.5, -6)) +
  scale_y_continuous(limits = c(4.5, 6.5), labels = NULL) +
  labs(x = NULL, y = NULL)

#Universal:
gene2439_12dpi_16 <- FeaturePlot(HC_2.srt, features = "gene2439", pt.size = 0.5, cols = c("grey", "grey")) +
  theme(axis.text.x = element_text(size = 13, angle = 45, vjust = 0.5)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(legend.position = 'none') +
  ggtitle(NULL) +
  scale_x_continuous(limits = c(-7.5, -6)) +
  scale_y_continuous(limits = c(4.5, 6.5), labels = NULL) +
  labs(x = NULL, y = NULL)

gene452_12dpi_16 <- FeaturePlot(HC_2.srt, features = "gene452", pt.size = 0.5) +
  theme(axis.text.x = element_text(size = 13, angle = 45, vjust = 0.5)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(legend.position = 'none') +
  ggtitle(NULL) +
  scale_x_continuous(limits = c(-7.5, -6)) +
  scale_y_continuous(limits = c(4.5, 6.5), labels = NULL) +
  labs(x = NULL, y = NULL)







#COMBINE ALL PLOTS

plot <- gene9384_4dpi_14 + gene6205_4dpi_14 + gene2439_4dpi_14 + gene452_4dpi_14 + gene9384_4dpi_16 + gene6205_4dpi_16 + gene2439_4dpi_16 + gene452_4dpi_16 + patchwork::plot_layout(ncol = 4)

setwd("/Users/emilyfitzmeyer/Desktop/")
ggsave("HC_plot.png", plot = plot, device = png(), scale = 1, width = 8, height = 8.5, dpi = 300)
dev.off()




plot2 <- clusters_4dpi + clusters_12dpi + patchwork::plot_layout(ncol = 1)
ggsave("clusters.png", plot = plot2, device = png(), scale = 1, width = 8.5, height = 10, dpi = 300)

