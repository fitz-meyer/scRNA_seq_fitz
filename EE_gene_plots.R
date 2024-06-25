library(tidyverse)
library(Seurat)
library(gridExtra)
library(ggplot2)
library(scCustomize)


merged.srt <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq_pub/mergeDPI_oldNorm/alldpi_oldNorm_nmn_res0.6.rds")
split.srt <- SplitObject(merged.srt, split.by = "dpi")
dpi.srt <- split.srt$"12dpi"
EE.srt <- subset(dpi.srt, idents = 12)

DimPlot(EE.srt, reduction = 'umap')
#EE_plot_breaks <- c( "Tachykinin", "IA2", "SCG5 (7B2 precursor)", "Syt1", "Syt4", "Syt6", "nSyb", "Syx1A" )
#EE_plot_features <- c("gene13952", "gene8024", "gene13213", "gene2866",  "gene3638", "gene8406", "gene7584", "gene5192")

#make binary version of "data" slot
bin_dat_EE <- GetAssayData(EE.srt[["RNA"]], slot = "data")
bin_dat_EE[bin_dat_EE>0] <- 1
EE.srt[["binary"]] <- CreateAssayObject(data = bin_dat_EE)
DefaultAssay(EE.srt) <- "binary"

#TOP ROW
gene13952 <- FeaturePlot(EE.srt, features = "gene13952", pt.size = 0.5) &
  theme(axis.text.x = element_text(size = 15)) &
  theme(axis.text.y = element_text(size = 15)) &
  theme(legend.position = 'none') &
  ggtitle("Tk (receptor)") &
  scale_x_discrete(labels = 'none') &
  labs(x = NULL, y = NULL)

gene8024 <- FeaturePlot(EE.srt, features = "gene8024", pt.size = 0.5) &
  theme(axis.text.x = element_text(size = 15)) &
  theme(axis.text.y = element_text(size = 15)) &
  theme(legend.position = 'none') &
  ggtitle("IA2") &
  scale_x_discrete(labels = 'none') &
  scale_y_discrete(labels = 'none') &
  labs(x = NULL, y = NULL)

gene13213 <- FeaturePlot(EE.srt, features = "gene13213", pt.size = 0.5) &
  theme(axis.text.x = element_text(size = 15)) &
  theme(axis.text.y = element_text(size = 15)) &
  theme(legend.position = 'none') &
  ggtitle("SCG5") &
  scale_x_discrete(labels = 'none') &
  scale_y_discrete(labels = 'none') &
  labs(x = NULL, y = NULL)

#MIDDLE ROW
gene2866 <- FeaturePlot(EE.srt, features = "gene2866", pt.size = 0.5) &
  theme(axis.text.x = element_text(size = 15)) &
  theme(axis.text.y = element_text(size = 15)) &
  theme(legend.position = 'none') &
  ggtitle("Syt1") &
  scale_x_discrete(labels = 'none') &
  labs(x = NULL, y = NULL)

gene3638 <- FeaturePlot(EE.srt, features = "gene3638", pt.size = 0.5) &
  theme(axis.text.x = element_text(size = 15)) &
  theme(axis.text.y = element_text(size = 15)) &
  theme(legend.position = 'none') &
  ggtitle("Syt4") &
  scale_x_discrete(labels = 'none') &
  scale_y_discrete(labels = 'none') &
  labs(x = NULL, y = NULL)

gene8406 <- FeaturePlot(EE.srt, features = "gene8406", pt.size = 0.5) &
  theme(axis.text.x = element_text(size = 15)) &
  theme(axis.text.y = element_text(size = 15)) &
  theme(legend.position = 'none') &
  ggtitle("Syt6") &
  scale_x_discrete(labels = 'none') &
  scale_y_discrete(labels = 'none') &
  labs(x = NULL, y = NULL)

#BOTTOM ROW
gene7584 <- FeaturePlot(EE.srt, features = "gene7584", pt.size = 0.5) &
  theme(axis.text.x = element_text(size = 15)) &
  theme(axis.text.y = element_text(size = 15)) &
  theme(legend.position = 'none') &
  ggtitle("nSyb") &
  labs(x = NULL, y = NULL)

gene5192 <- FeaturePlot(EE.srt, features = "gene5192", pt.size = 0.5) &
  theme(axis.text.x = element_text(size = 15)) &
  theme(axis.text.y = element_text(size = 15)) &
  theme(legend.position = 'none') &
  ggtitle("Syx1A") &
  scale_y_discrete(labels = 'none') &
  labs(x = NULL, y = NULL)

gene7998 <- FeaturePlot(EE.srt, features = "gene7998", pt.size = 0.5) &
  theme(axis.text.x = element_text(size = 15)) &
  theme(axis.text.y = element_text(size = 15)) &
  theme(legend.position = 'none') &
  ggtitle("NEUROD6") &
  scale_y_discrete(labels = 'none') &
  labs(x = NULL, y = NULL)

# #TEST
# gene1511 <- FeaturePlot(EE.srt, features = "gene1511", pt.size = 0.5) &
#   theme(axis.text.x = element_text(size = 15)) &
#   theme(axis.text.y = element_text(size = 15)) &
#   theme(legend.position = 'none') &
#   ggtitle("Burs") &
#   scale_y_discrete(labels = 'none') &
#   labs(x = NULL, y = NULL)
# 
# gene12241 <- FeaturePlot(EE.srt, features = "gene12241", pt.size = 0.5) &
#   theme(axis.text.x = element_text(size = 15)) &
#   theme(axis.text.y = element_text(size = 15)) &
#   theme(legend.position = 'none') &
#   ggtitle("ITP") &
#   scale_y_discrete(labels = 'none') &
#   labs(x = NULL, y = NULL)
# 
# gene2206 <- FeaturePlot(EE.srt, features = "gene2206", pt.size = 0.5) &
#   theme(axis.text.x = element_text(size = 15)) &
#   theme(axis.text.y = element_text(size = 15)) &
#   theme(legend.position = 'none') &
#   ggtitle("sNPF") &
#   scale_y_discrete(labels = 'none') &
#   labs(x = NULL, y = NULL)

plot <- gene13952 + gene8024 + gene13213 + gene2866 + gene3638 + gene8406 + gene7584 + gene5192 + gene7998 + patchwork::plot_layout(ncol = 3)

setwd("/Users/emilyfitzmeyer/Desktop/")
ggsave("EE12_neurpep_plots.png", plot = plot, device = png(), scale = 1, width = 7, height = 7, dpi = 300)
dev.off()

