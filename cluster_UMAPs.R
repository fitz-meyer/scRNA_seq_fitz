library(tidyverse)
library(Seurat)
library(gridExtra)
library(ggplot2)
library(scales)

idents <- levels(merged.srt@active.ident)
my_color_palette <- hue_pal()(length(idents))

#4dpi
merged.srt <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq_pub/mergeDPI_oldNorm/alldpi_oldNorm_nmn_res0.6.rds")
new_idents <- c("EC-like-1", "1", "EC-like-2", "ISC/EB", "4", "EC", "EC-like-3", "VM-1", "8", "VM-2", "ISC/EB-prol", "cardia-1", "EE", "cardia-2", "HC-1", "MT", "HC-2", "17", "cardia-prol", "19")
merged.srt <- Rename_Clusters(merged.srt, new_idents = new_idents)

split.srt <- SplitObject(merged.srt, split.by = "dpi")
srt_4dpi <- split.srt$"4dpi"
srt_12dpi <- split.srt$"12dpi"
srt_4dpi <- SplitObject(srt_4dpi, split.by = "condition")
srt_12dpi <- SplitObject(srt_12dpi, split.by = "condition")

clusters_4dpi_mock <- DimPlot(srt_4dpi$mock, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.2, label.size = 4.5, 
                              cols = c("EC-like-1" = "#F8766D", "1" = "#EA8331", "EC-like-2" = "#D89000", "ISC/EB" = "#C09B00", 
                                       "4" = "#A3A500", "EC" = "#7CAE00", "EC-like-3" = "#39B600", "VM-1" = "#00BB4E", "8" = "#00BF7D",
                                       "VM-2" = "#00C1A3", "ISC/EB-prol" = "#00BFC4", "cardia-1" = "#00BAE0", "EE" = "#00B0F6", 
                                       "cardia-2" = "#9590FF", "HC-1" = "#E76BF3", "MT" = "#C77CFF", "HC-2" = "#35A2FF", "HC-2" = "#FA62DB", 
                                       "17" = "#FF62BC", "cardia-prol" = "#FF6A98", "19" = "darkgray")) +
  theme(axis.text.x = element_text(size = 17)) +
  theme(axis.text.y = element_text(size = 17)) +
  theme(legend.position = 'none') +
  labs(x = NULL, y = "4dpi") +
  theme(axis.title = element_text(size = 17))

clusters_4dpi_WNV <- DimPlot(srt_4dpi$wnv, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.2, label.size = 4.5, 
                             cols = c("EC-like-1" = "#F8766D", "1" = "#EA8331", "EC-like-2" = "#D89000", "ISC/EB" = "#C09B00", 
                                      "4" = "#A3A500", "EC" = "#7CAE00", "EC-like-3" = "#39B600", "VM-1" = "#00BB4E", "8" = "#00BF7D",
                                      "VM-2" = "#00C1A3", "ISC/EB-prol" = "#00BFC4", "cardia-1" = "#00BAE0", "EE" = "#00B0F6", 
                                      "cardia-2" = "#9590FF", "HC-1" = "#E76BF3", "MT" = "#C77CFF", "HC-2" = "#35A2FF", "HC-2" = "#FA62DB", 
                                      "17" = "#FF62BC", "cardia-prol" = "#FF6A98", "19" = "darkgray")) +
  theme(axis.text.x = element_text(size = 17)) +
  theme(axis.text.y = element_text(size = 17)) +
  theme(legend.position = 'none') +
  labs(x = NULL, y = NULL) +
  scale_y_continuous(labels = NULL)

clusters_4dpi_merge <- DimPlot(split.srt$"4dpi", reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.2, label.size = 4.5,
                               cols = c("EC-like-1" = "#F8766D", "1" = "#EA8331", "EC-like-2" = "#D89000", "ISC/EB" = "#C09B00", 
                                        "4" = "#A3A500", "EC" = "#7CAE00", "EC-like-3" = "#39B600", "VM-1" = "#00BB4E", "8" = "#00BF7D",
                                        "VM-2" = "#00C1A3", "ISC/EB-prol" = "#00BFC4", "cardia-1" = "#00BAE0", "EE" = "#00B0F6", 
                                        "cardia-2" = "#9590FF", "HC-1" = "#E76BF3", "MT" = "#C77CFF", "HC-2" = "#35A2FF", "HC-2" = "#FA62DB", 
                                        "17" = "#FF62BC", "cardia-prol" = "#FF6A98", "19" = "darkgray")) +
  theme(axis.text.x = element_text(size = 17)) +
  theme(axis.text.y = element_text(size = 17)) +
  theme(legend.position = 'none') +
  scale_y_continuous(labels = NULL) +
  labs(x = NULL, y = NULL)

plot_1 <- clusters_4dpi_mock + clusters_4dpi_WNV + clusters_4dpi_merge + patchwork::plot_layout(ncol = 3)

#12dpi

clusters_12dpi_mock <- DimPlot(srt_12dpi$mock, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.2, label.size = 4.5,
                               cols = c("EC-like-1" = "#F8766D", "1" = "#EA8331", "EC-like-2" = "#D89000", "ISC/EB" = "#C09B00", 
                                        "4" = "#A3A500", "EC" = "#7CAE00", "EC-like-3" = "#39B600", "VM-1" = "#00BB4E", "8" = "#00BF7D",
                                        "VM-2" = "#00C1A3", "ISC/EB-prol" = "#00BFC4", "cardia-1" = "#00BAE0", "EE" = "#00B0F6", 
                                        "cardia-2" = "#9590FF", "HC-1" = "#E76BF3", "MT" = "#C77CFF", "HC-2" = "#35A2FF", "HC-2" = "#FA62DB", 
                                        "17" = "#FF62BC", "cardia-prol" = "#FF6A98", "19" = "darkgray")) +
  theme(axis.text.x = element_text(size = 17)) +
  theme(axis.text.y = element_text(size = 17)) +
  theme(legend.position = 'none') +
  labs(x = NULL, y = "12dpi") +
  theme(axis.title = element_text(size = 17))

clusters_12dpi_WNV <- DimPlot(srt_12dpi$wnv, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.2, label.size = 4.5,
                              cols = c("EC-like-1" = "#F8766D", "1" = "#EA8331", "EC-like-2" = "#D89000", "ISC/EB" = "#C09B00", 
                                       "4" = "#A3A500", "EC" = "#7CAE00", "EC-like-3" = "#39B600", "VM-1" = "#00BB4E", "8" = "#00BF7D",
                                       "VM-2" = "#00C1A3", "ISC/EB-prol" = "#00BFC4", "cardia-1" = "#00BAE0", "EE" = "#00B0F6", 
                                       "cardia-2" = "#9590FF", "HC-1" = "#E76BF3", "MT" = "#C77CFF", "HC-2" = "#35A2FF", "HC-2" = "#FA62DB", 
                                       "17" = "#FF62BC", "cardia-prol" = "#FF6A98", "19" = "darkgray")) +
  theme(axis.text.x = element_text(size = 17)) +
  theme(axis.text.y = element_text(size = 17)) +
  theme(legend.position = 'none') +
  scale_y_continuous(labels = NULL) +
  labs(x = NULL, y = NULL)

clusters_12dpi_merge <- DimPlot(split.srt$"12dpi", reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.2, label.size = 4.5,
                                cols = c("EC-like-1" = "#F8766D", "1" = "#EA8331", "EC-like-2" = "#D89000", "ISC/EB" = "#C09B00", 
                                         "4" = "#A3A500", "EC" = "#7CAE00", "EC-like-3" = "#39B600", "VM-1" = "#00BB4E", "8" = "#00BF7D",
                                         "VM-2" = "#00C1A3", "ISC/EB-prol" = "#00BFC4", "cardia-1" = "#00BAE0", "EE" = "#00B0F6", 
                                         "cardia-2" = "#9590FF", "HC-1" = "#E76BF3", "MT" = "#C77CFF", "HC-2" = "#35A2FF", "HC-2" = "#FA62DB", 
                                         "17" = "#FF62BC", "cardia-prol" = "#FF6A98", "19" = "darkgray")) +
  theme(axis.text.x = element_text(size = 17)) +
  theme(axis.text.y = element_text(size = 17)) +
  theme(legend.position = 'none') +
  scale_y_continuous(labels = NULL) +
  labs(x = NULL, y = NULL)

plot_2 <- clusters_12dpi_mock + clusters_12dpi_WNV + clusters_12dpi_merge + patchwork::plot_layout(ncol = 3)

plot_3 <- clusters_4dpi_mock + clusters_4dpi_WNV + clusters_4dpi_merge + clusters_12dpi_mock + clusters_12dpi_WNV + clusters_12dpi_merge + patchwork::plot_layout(ncol = 3)

ggsave("clusters_ALL.png", plot = plot_3, device = png(), scale = 1, width = 12.5, height = 8.5, dpi = 300)
dev.off()


