#DotPlots

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

split.srt <- SplitObject(merged.srt, split.by = "dpi")
dpi4.srt <- split.srt$"4dpi"
dpi12.srt <- split.srt$"12dpi"
dpi4.srt <- SplitObject(dpi4.srt, split.by = "condition")
dpi12.srt <- SplitObject(dpi12.srt, split.by = "condition")



#CLUSTER MARKER DOTPLOT
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
#levels(merged.srt@meta.data$RNA_snn_res.0.6)
# 
# Idents(merged.srt) <- "cell_type"
# Idents(merged.srt) <- factor(merged.srt@active.ident, sort(levels(merged.srt@active.ident)))

#reorder y axis to list cell types alphabetically/group types together (no number weirdness)
my_levels <- c("VM-2", "VM-1", "MT", "ISC/EB-prol", "ISC/EB", "HC-2", "HC-1", "EE", "EC-like-3", "EC-like-2", "EC-like-1", "EC", "cardia-prol", "cardia-2", "cardia-1", "19", "17", "8", "4", "1")
factor(Idents(merged.srt), levels = my_levels)
Idents(merged.srt) <- factor(Idents(merged.srt), levels = my_levels)

DotPlot_scCustom(merged.srt, features = dotplot_genes) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = dotplot_break_labels)

setwd("/Users/emilyfitzmeyer/Desktop/")
ggsave("cl_marker_dotplot_alldpi.png", plot = last_plot(), device = png(), scale = 1, width = 9, height = 6.5, dpi = 300)
dev.off()

#cl4 (high vRNA cluster) dotplot
dotplot_genes <- c("nbisL1-trna-9", "MT-nbis-gene-3", "MT-nbis-gene-2", "nbisL1-trna-7", "nbisL1-trna-10",
                   "nbisL1-trna-6", "nbisL1-trna-16", "nbisL1-trna-3", "nbisL1-trna-17", "MT-nbis-gene-5",
                   "nbisL1-trna-20", "nbisL1-trna-19", "nbis-gene-2-utr")

dotplot_break_labels <- c("MT-tRNA-Asp", "COX3", "ATP6", "MT-tRNA-Leu", "MT-tRNA-Gly",
                          "MT-tRNA-Tyr", "MT-tRNA-Phe", "MT-tRNA-Met", "MT-tRNA-His", "CYTB",
                          "MT-tRNA-Ser", "MT-tRNA-Pro", "WNV 5' UTR")

#reorder y axis to list cell types alphabetically/group types together (no number weirdness)
my_levels <- c("VM-2", "VM-1", "MT", "ISC/EB-prol", "ISC/EB", "HC-2", "HC-1", "EE", "EC-like-3", "EC-like-2", "EC-like-1", "EC", "cardia-prol", "cardia-2", "cardia-1", "19", "17", "8", "4", "1")
factor(Idents(dpi4.srt$wnv), levels = my_levels)
Idents(dpi4.srt$wnv) <- factor(Idents(dpi4.srt$wnv), levels = my_levels)

dpi4_cl4_dotplot <- DotPlot_scCustom(dpi4.srt$wnv, features = dotplot_genes) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = 'none') +
  scale_x_discrete(labels = dotplot_break_labels)

#reorder y axis to list cell types alphabetically/group types together (no number weirdness)
my_levels <- c("VM-2", "VM-1", "MT", "ISC/EB-prol", "ISC/EB", "HC-2", "HC-1", "EE", "EC-like-3", "EC-like-2", "EC-like-1", "EC", "cardia-prol", "cardia-2", "cardia-1", "19", "17", "8", "4", "1")
factor(Idents(dpi12.srt$wnv), levels = my_levels)
Idents(dpi12.srt$wnv) <- factor(Idents(dpi12.srt$wnv), levels = my_levels)

dpi12_cl4_dotplot <- DotPlot_scCustom(dpi12.srt$wnv, features = dotplot_genes) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = dotplot_break_labels)

plot <- dpi4_cl4_dotplot + dpi12_cl4_dotplot + patchwork::plot_layout(ncol = 2)

ggsave("cl4_marker_dotplot.png", plot = plot, device = png(), scale = 1, width = 10.5, height = 5, dpi = 300)
dev.off()

#ISC/EB clusters dotplot
dotplot_genes <- c("gene6060", "gene4986", "gene10089",	"gene6115",	"gene6292",	"gene1252",	"gene6564",
                   "gene2793", "gene6552", "gene5534", "gene4419")

dotplot_break_labels <- c("DUF4803", "gene4986", "cecropin", "cecropin", "CUTA", "GSTD11", "awd",
                          "NHP2L1", "histone H2A", "RAN", "BIRC5")

DotPlot_scCustom(merged.srt, features = dotplot_genes) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = dotplot_break_labels)

ggsave("ISC_marker_dotplot_4dpi.png", plot = last_plot(), device = png(), scale = 1, width = 6, height = 4, dpi = 300)
dev.off()





