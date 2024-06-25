#UMAPs and FeaturePlots

library(singleCellTK)
library(tidyverse)
library(Seurat)
library(scCustomize)
library(gridExtra)
library(EnhancedVolcano)
library(clustree)

#Read in merged.srt
merged.srt <- readRDS("/Users/emilyfitzmeyer/Desktop/scRNAseq_pub/mergeDPI_oldNorm/alldpi_oldNorm_nmn_res0.6.rds")
#rename clusters
new_idents <- c("EC-like-1", "1", "EC-like-2", "ISC/EB", "4", "EC", "EC-like-3", "VM-1", "8", "VM-2", "ISC/EB-prol", "cardia-1", "EE", "cardia-2", "HC-1", "MT", "HC-2", "17", "cardia-prol", "19")
merged.srt <- Rename_Clusters(merged.srt, new_idents = new_idents)
#add cell type info to metadata
ident_num <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19")
cell_type <- c("EC-like-1", "1", "EC-like-2", "ISC/EB", "4", "EC", "EC-like-3", "VM-1", "8", "VM-2", "ISC/EB-prol", "cardia-1", "EE", "cardia-2", "HC-1", "MT", "HC-2", "17", "cardia-prol", "19")
tib <- tibble(ident_num, cell_type)
merged.srt <- Add_Sample_Meta(merged.srt, meta_data = tib, join_by_seurat = "seurat_clusters", join_by_meta = "ident_num")

my_levels <- c("VM-2", "VM-1", "MT", "ISC/EB-prol", "ISC/EB", "HC-2", "HC-1", "EE", "EC-like-3", "EC-like-2", "EC-like-1", "EC", "cardia-prol", "cardia-2", "cardia-1", "19", "17", "8", "4", "1")
factor(Idents(merged.srt), levels = my_levels)
Idents(merged.srt) <- factor(Idents(merged.srt), levels = my_levels)

split.srt <- SplitObject(merged.srt, split.by = "dpi")
srt_4dpi <- split.srt$"4dpi"
srt_12dpi<- split.srt$"12dpi"
srt_4dpi <- SplitObject(srt_4dpi, split.by = "condition")
srt_12dpi <- SplitObject(srt_12dpi, split.by = "condition")

#show sample and condition induced impact on cell grouping:
sample_plot_4dpi <- DimPlot(srt_4dpi, reduction = 'umap', group.by = 'sample', pt.size = 0.1) &
  theme(axis.text.x = element_text(size = 15)) &
  theme(axis.text.y = element_text(size = 15)) &
  scale_x_continuous(labels = NULL) &
  ggtitle("Sample Grouping") &
  labs(x = NULL)

sample_plot_12dpi <- DimPlot(srt_12dpi, reduction = 'umap', group.by = 'sample', pt.size = 0.1) &
  theme(axis.text.x = element_text(size = 15)) &
  theme(axis.text.y = element_text(size = 15)) &
  ggtitle(NULL)

condition_plot_4dpi <- DimPlot(srt_4dpi, reduction = 'umap', group.by = 'condition',
                          cols = c('deepskyblue', 'salmon'), pt.size = 0.1) &
  theme(axis.text.x = element_text(size = 15)) &
  theme(axis.text.y = element_text(size = 15)) &
  scale_x_continuous(labels = NULL) &
  scale_y_continuous(labels = NULL) &
  ggtitle("Condition Grouping") &
  labs(x = NULL, y = NULL)

condition_plot_12dpi <- DimPlot(srt_12dpi, reduction = 'umap', group.by = 'condition',
                               cols = c('deepskyblue', 'salmon'), pt.size = 0.1) &
  theme(axis.text.x = element_text(size = 15)) &
  theme(axis.text.y = element_text(size = 15)) &
  scale_y_continuous(labels = NULL) &
  ggtitle(NULL) &
  labs(y = NULL)

plot <- sample_plot_4dpi + condition_plot_4dpi + sample_plot_12dpi + condition_plot_12dpi + patchwork::plot_layout(ncol = 2)

setwd("/Users/emilyfitzmeyer/Desktop/")
ggsave("batch_plots.png", plot = plot, device = png(), scale = 1, width = 8.5, height = 6.5, dpi = 300)
dev.off()


#plot WNV 5' UTR:
WNV.utr.4dpi.mock <- FeaturePlot(srt_4dpi$mock, features = "nbis-gene-2-utr", pt.size = 0.5, cols = c("grey", "grey")) +
  theme(axis.text.x = element_text(size = 17)) +
  theme(axis.text.y = element_text(size = 17)) +
  theme(legend.position = 'none') +
  scale_x_discrete(labels = 'none') +
  ggtitle("Mock") +
  labs(x = NULL) +
  theme(axis.title.y = element_text(vjust = -2.5))

WNV.utr.4dpi.wnv <- FeaturePlot(srt_4dpi$wnv, features = "nbis-gene-2-utr", pt.size = 0.5) +
  theme(axis.text.x = element_text(size = 17)) +
  theme(axis.text.y = element_text(size = 17)) +
  theme(legend.position = c(0.92, 0.92)) +
  ggtitle("WNV") +
  scale_x_discrete(labels = 'none') +
  scale_y_discrete(labels = 'none') +
  labs(x = NULL, y = NULL)

WNV.utr.12dpi.mock <- FeaturePlot(srt_12dpi$mock, features = "nbis-gene-2-utr", pt.size = 0.5, cols = c("grey", "grey")) +
  theme(axis.text.x = element_text(size = 17)) +
  theme(axis.text.y = element_text(size = 17)) +
  theme(legend.position = 'none') +
  ggtitle(NULL) +
  theme(axis.title.y = element_text(vjust = -2.5))

WNV.utr.12dpi.wnv <- FeaturePlot(srt_12dpi$wnv, features = "nbis-gene-2-utr", pt.size = 0.5) +
  theme(axis.text.x = element_text(size = 17)) +
  theme(axis.text.y = element_text(size = 17)) +
  theme(legend.position = 'none') +
  ggtitle(NULL) +
  scale_y_discrete(labels = 'none') +
  labs(y = NULL)

plot <- WNV.utr.4dpi.mock + WNV.utr.4dpi.wnv + WNV.utr.12dpi.mock + WNV.utr.12dpi.wnv + patchwork::plot_layout(ncol = 2)

setwd("/Users/emilyfitzmeyer/Desktop/")
ggsave("utr_load_featureplots.png", plot = plot, device = png(), scale = 1, width = 7.5, height = 7.5, dpi = 300)
dev.off()


#view cell cycle features__________________________________________________________________________________
#S phase
FeaturePlot(merged.srt, features = c("gene7803", "gene1526", "gene2166", "gene577", "gene9562"))
ggsave("sPhase_bothdpi_plots.png", plot = last_plot(), device = png(), scale = 1, width = 8, height = 10, dpi = 300)
dev.off()

#G2/M phase
FeaturePlot(merged.srt, features = c("gene13155", "gene9499", "gene7697", "gene4233", "gene11073", "gene13413", "gene13452", "gene5102", "gene9506", "gene6438", "gene7563", "gene9947"))
ggsave("g2mPhase_bothdpi_plots.png", plot = last_plot(), device = png(), scale = 1, width = 15, height = 10, dpi = 300)
dev.off()

#view apoptotic genes__________________________________________________________________________________
FeaturePlot(merged.srt, features = c("gene2212", "gene13832", "gene9932", "gene11674", "gene4419", "gene9707", "gene1473", "gene13123", "gene14035"))
ggsave("apoptotic_bothdpi_plots.png", plot = last_plot(), device = png(), scale = 1, width = 12, height = 10, dpi = 300)
dev.off()

#ISC plots
klu_Fplot <- FeaturePlot(merged.srt, features = "gene13104", pt.size = 0.2) +
  theme(axis.text.x = element_text(size = 17)) +
  theme(axis.text.y = element_text(size = 17)) +
  theme(legend.position = c(0.1, 0.95)) +
  # scale_x_continuous(limits = c(-6, 4)) +
  # scale_y_continuous(limits = c(-14, -7)) +
  theme(legend.text = element_text(size = 17)) +
  ggtitle("Klu")

klu_Vplot <- VlnPlot(merged.srt, features = "gene13104") +
  theme(legend.position = 'none') + 
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 17)) +
  labs(x = NULL) +
  ggtitle("Klu")

plot <- klu_Fplot + klu_Vplot + patchwork::plot_layout(ncol = 2)

ggsave("klu_bothdpi_plots.png", plot = plot, device = png(), scale = 1, width = 9, height = 5, dpi = 300)
dev.off()

#ISC-prol plots
pcna_Fplot <- FeaturePlot(merged.srt, features = "gene2167", pt.size = 0.5) +
  theme(axis.text.x = element_text(size = 17)) +
  theme(axis.text.y = element_text(size = 17)) +
  scale_x_continuous(limits = c(-6, 4), labels = NULL) +
  scale_y_continuous(limits = c(-14, -7)) +
  theme(legend.position = c(0.9, 0.23)) +
  theme(legend.text = element_text(size = 17)) +
  theme(legend.position = c(0.85, 0.95)) +
  labs(x = NULL) +
  ggtitle("PCNA")

pcna_Vplot <- VlnPlot(merged.srt, features = "gene2167") +
  theme(legend.position = 'none') + 
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 17)) +
  scale_x_discrete(labels = NULL) +
  labs(x = NULL) +
  ggtitle("PCNA")

aurka_Fplot <- FeaturePlot(merged.srt, features = "gene13413", pt.size = 0.5) +
  theme(axis.text.x = element_text(size = 17)) +
  theme(axis.text.y = element_text(size = 17)) +
  scale_x_continuous(limits = c(-6, 4), labels = NULL) +
  scale_y_continuous(limits = c(-14, -7)) +
  theme(legend.position = c(0.9, 0.23)) +
  theme(legend.text = element_text(size = 17)) +
  theme(legend.position = 'none') +
  labs(x = NULL) +
  ggtitle("AURKA")

aurka_Vplot <- VlnPlot(merged.srt, features = "gene13413") +
  theme(legend.position = 'none') + 
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 17)) +
  scale_x_discrete(labels = NULL) +
  labs(x = NULL) +
  ggtitle("AURKA")

aurkb_Fplot <- FeaturePlot(merged.srt, features = "gene12926", pt.size = 0.5) +
  theme(axis.text.x = element_text(size = 17)) +
  theme(axis.text.y = element_text(size = 17)) +
  scale_x_continuous(limits = c(-6, 4)) +
  scale_y_continuous(limits = c(-14, -7)) +
  theme(legend.position = c(0.9, 0.23)) +
  theme(legend.text = element_text(size = 17)) +
  theme(legend.position = 'none') +
  ggtitle("AURKB")

aurkb_Vplot <- VlnPlot(merged.srt, features = "gene12926") +
  theme(legend.position = 'none') + 
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 17)) +
  labs(x = NULL) +
  ggtitle("AURKB")

plot <- pcna_Fplot + pcna_Vplot + aurka_Fplot + aurka_Vplot + aurkb_Fplot + aurkb_Vplot + patchwork::plot_layout(ncol = 2)

ggsave("prol_marker_bothdpi_plots.png", plot = plot, device = png(), scale = 1, width = 9, height = 11, dpi = 300)
dev.off()


#PROX1 plots
prox1_Fplot <- FeaturePlot(merged.srt, features = "gene10632") +
  theme(axis.text.x = element_text(size = 17)) +
  theme(axis.text.y = element_text(size = 17)) +
  theme(legend.position = c(0.1, 0.95)) +
  ggtitle("PROX1")

prox1_Vplot <- VlnPlot(merged.srt, features = "gene10632") + 
  theme(legend.position = 'none') + 
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 17)) +
  labs(x = NULL) +
  ggtitle("PROX1")

plot <- prox1_Fplot + prox1_Vplot + patchwork::plot_layout(ncol = 2)
  
ggsave("PROX1_plots_bothdpi.png", plot = last_plot(), device = png(), scale = 1, width = 9, height = 5, dpi = 300)
dev.off()


#VM plots
actb_Fplot <- FeaturePlot(merged.srt, features = "gene2985") +
  theme(axis.text.x = element_text(size = 17)) +
  theme(axis.text.y = element_text(size = 17)) +
  scale_x_continuous(labels = NULL) +
  theme(legend.position = c(0.05, 0.95)) +
  labs(x = NULL) +
  ggtitle("ACTB")

mhc_Fplot <- FeaturePlot(merged.srt, features = "gene11104") +
  theme(axis.text.x = element_text(size = 17)) +
  theme(axis.text.y = element_text(size = 17)) +
  scale_x_continuous(labels = NULL) +
  theme(legend.position = 'none') +
  labs(x = NULL) +
  ggtitle("Mhc/Mhc1")

mlc_Fplot <- FeaturePlot(merged.srt, features = "gene9804") +
  theme(axis.text.x = element_text(size = 17)) +
  theme(axis.text.y = element_text(size = 17)) +
  theme(legend.position = 'none') +
  ggtitle("Mlc2")


actb_Vplot <- VlnPlot(merged.srt, features = "gene2985") + 
  theme(legend.position = 'none') + 
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 17)) +
  scale_x_discrete(labels = NULL) +
  labs(x = NULL) +
  ggtitle("ACTB")

mhc_Vplot <- VlnPlot(merged.srt, features = "gene11104") + 
  theme(legend.position = 'none') + 
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 17)) +
  scale_x_discrete(labels = NULL) +
  labs(x = NULL) +
  ggtitle("Mhc/Mhc1")

mlc2_Vplot <- VlnPlot(merged.srt, features = "gene9804") + 
  theme(legend.position = 'none') + 
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 17)) +
  labs(x = NULL) +
  ggtitle("Mlc2")


plot <- actb_Fplot + actb_Vplot + mhc_Fplot + mhc_Vplot + mlc_Fplot + mlc2_Vplot + patchwork::plot_layout(ncol = 2)

ggsave("VM_markers_bothdpi.png", plot = plot, device = png(), scale = 1, width = 9, height = 11, dpi = 300)
dev.off()


#EC plots
FeaturePlot(merged.srt, features = "gene11957")




