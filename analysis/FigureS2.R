library(Seurat)
library(tidyverse)
library(ggpubr)
library(MetaNeighbor)
library(clustree)
library(SingleCellExperiment)
library(ComplexHeatmap)

young_obj <- readRDS("/Volumes/LACIE/RPE_scRNA/PUBLISHED/H9_RPE_YOUNG_Analysed_Object.rds")
aged_obj <- readRDS("/Volumes/LACIE/RPE_scRNA/PUBLISHED/H9_RPE_AGED_Analysed_Object.rds")
integrated_obj <- readRDS("/Volumes/LACIE/RPE_scRNA/PUBLISHED/H9_RPE_IntegratedObject.rds")

young_clustree <- clustree(young_obj, prefix = "SCT_snn_res.")
aged_clustree <- clustree(aged_obj, prefix = "SCT_snn_res.")
integrated_clustree <- clustree(integrated_obj, prefix = "integrated_snn_res.")

clustree_plot <- ggarrange(young_clustree, aged_clustree, ncol = 2)

library(MetaNeighbor)
library(SingleCellExperiment)

# Prepare data for conversion to SingleCellExperiment
metadata <- FetchData(integrated_obj, vars = c("orig.ident", "condition", "seurat_clusters"))
metadata$cluster <- paste0(metadata$condition, "_", metadata$seurat_clusters)
integrated_obj[["cluster"]] <- metadata$cluster
sce_obj <- as.SingleCellExperiment(integrated_obj)

# Run unsupervised metaneighbor based on top variable genes
variable_genes <- variableGenes(sce_obj, exp_labels = sce_obj$condition)
celltype_nv = MetaNeighborUS(var_genes = variable_genes,
                             dat = sce_obj,
                             study_id = sce_obj$condition,
                             cell_type = sce_obj$seurat_clusters, fast_version = TRUE)

colnames(celltype_nv) <- gsub("^Control", "Young", colnames(celltype_nv))
rownames(celltype_nv) <- gsub("^Control", "Young", rownames(celltype_nv))

cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
breaks = seq(0, 1, length=101)

auroc_matrix <- celltype_nv[paste0("Aged|", as.character(0:11)), 
                            paste0("Young|", as.character(0:11))]
auroc_plot <- Heatmap(auroc_matrix, col = cols, name = "AUROC",
                      rect_gp = gpar(col = "white", lwd = 2))

ggsave("FigureS2A.pdf", clustree_plot, width = 220, height = 300, units = "mm", dpi = 600)

pdf("FigureS2B.pdf", width = 7, height = 5)
draw(auroc_plot)
dev.off()

ggsave("FigureS2Bii.pdf", integrated_clustree, width = 180, height = 280, units = "mm", dpi = 600)

dev.off()

