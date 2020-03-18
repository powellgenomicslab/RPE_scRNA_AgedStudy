library(Seurat)
library(tidyverse)
library(patchwork)
library(slingshot)
library(monocle3)
library(ggpubr)

seurat_obj <- readRDS("/Volumes/LACIE/RPE_scRNA/PUBLISHED/H9_RPE_IntegratedObject.rds")
monocle_obj <- readRDS("/Volumes/LACIE/RPE_scRNA/PUBLISHED/H9_RPE_Monocle3_Object.rds")
slingshot_obj <- readRDS("/Volumes/LACIE/RPE_scRNA/PUBLISHED/H9_AGED_RPE_SlingshotObject.rds")

seurat_obj$condition[seurat_obj$condition == "Control"] <- "Young"
saveRDS(seurat_obj, "/Volumes/LACIE/RPE_scRNA/PUBLISHED/H9_RPE_IntegratedObject.rds")

colour_df <- read_csv("~/Documents/Experiments/RPE_scRNA/Analysis/RevisedAnalysis/cluster_colours.csv")
colour_df$Cluster <- factor(colour_df$Cluster, levels = c(paste0("Young_", as.character(0:5)),
                                                          paste0("Aged_", as.character(0:5)),
                                                          paste0("Common_", as.character(1:6)))
)

colour_vec <- colour_df$Colour
names(colour_vec) <- colour_df$Cluster
colData(monocle_obj)$seurat_clusters <- seurat_obj[["seurat_clusters"]]$seurat_clusters

# Basic cluster characterisation
split_cluster_umap <- DimPlot(seurat_obj, reduction = "umap", group.by = "cluster",
                        cols = colour_vec, split.by = "condition") 

ggsave("ClusterByTimepoint_UMAP.pdf", split_cluster_umap, width = 11, height = 7, units = "in")
ggsave("ClusterByTimepoint_UMAP.tiff", split_cluster_umap, width = 11, height = 6.5, units = "in", dpi = 600)

cluster_umap <- DimPlot(seurat_obj, reduction = "umap", group.by = "cluster",
                              cols = colour_vec) 

ggsave("Cluster_UMAP.pdf", cluster_umap, width = 11, height = 7, units = "in")
ggsave("Cluster_UMAP.tiff", cluster_umap, width = 11, height = 6.5, units = "in", dpi = 600)


cluster_umap <- DimPlot(seurat_obj, reduction = "umap", group.by = "cluster",
                              cols = colour_vec) + ggtitle("Clusters") 

ggsave("Cluster_UMAP.pdf", cluster_umap, width = 6, height = 5, units = "in")
ggsave("Cluster_UMAP.tiff", cluster_umap, width = 6, height = 5, units = "in", dpi = 600)

timepoint_umap <- DimPlot(seurat_obj, reduction = "umap", group.by = "condition",
                          cols = list("Young" = "#e71837", "Aged" = "#00674f")) +
  ggtitle("Timepoint")
ggsave("Timepoint_UMAP.pdf", timepoint_umap, width = 6, height = 5, units = "in")
ggsave("Timepoint_UMAP.tiff", timepoint_umap, width = 6, height = 5, units = "in", dpi = 600)


# Trajectory analysis
colnames(colData(monocle_obj))[ncol(colData(monocle_obj))] <- "Cluster"
cluster_trajectory_plot <- plot_cells(monocle_obj, color_cells_by = "Cluster", 
                               reduction_method = "UMAP", 
                               show_trajectory_graph = TRUE,
                               label_cell_groups = FALSE,
                               label_leaves = FALSE,
                               label_branch_points = FALSE) +
  scale_color_manual(name = "Cluster", values = colour_vec) +
  ggtitle("Clusters", subtitle = "Monocle 3 Trajectory Analysis")

ggsave("ClusterTrajectory_UMAP.pdf", cluster_trajectory_plot, width = 6.5, height = 5.5, units = "in")
ggsave("ClusterTrajectory_UMAP.tiff", cluster_trajectory_plot, width = 6.5, height = 5.5, units = "in", dpi = 600)



colnames(colData(monocle_obj))[1] <- "Timepoint"
colData(monocle_obj)$Timepoint[colData(monocle_obj)$Timepoint == "Control"] <- "Young"
timepoint_trajectory_plot <- plot_cells(monocle_obj, color_cells_by = "Timepoint", 
                               reduction_method = "UMAP", 
                               show_trajectory_graph = FALSE,
                               label_cell_groups = FALSE) +
  scale_color_manual(name = "Timepoint", values = list("Young" = "#e71837", "Aged" = "#00674f")) +
  ggtitle("Clusters", subtitle = "Monocle 3 Trajectory Analysis")

ggsave("TimepointTrajectory_UMAP.pdf", timepoint_trajectory_plot, width = 6.5, height = 5.5, units = "in")
ggsave("TimepointTrajectory_UMAP.tiff", timepoint_trajectory_plot, width = 6.5, height = 5.5, units = "in", dpi = 600)

pseudotime_plot <- plot_cells(monocle_obj, color_cells_by = "pseudotime", 
                              show_trajectory_graph = FALSE,
                              label_cell_groups = FALSE, label_leaves = FALSE, 
                              label_branch_points = FALSE) +
  ggtitle("Pseudotime", subtitle = "Monocle 3 Trajectory Analysis")
ggsave("Pseudotime_Monocle3_UMAP.pdf", pseudotime_plot, width = 6.5, height = 5.5, units = "in")
ggsave("Pseudotime_Monocle3_UMAP.tiff", pseudotime_plot, width = 6.5, height = 5.5, units = "in", dpi = 600)


monocle3_data <- colData(monocle_obj)
monocle3_data$partition <- partitions(monocle_obj)
monocle3_data$pseudotime <- pseudotime(monocle_obj)
monocle3_data <- as.data.frame(monocle3_data)

monocle3_cluster_count <- as.data.frame(monocle3_data) %>% group_by(cluster, condition) %>% 
  summarize(n_cells = n_distinct(cell_barcode), partitions = n_distinct(partition))

pseudotime_by_cluster <- ggplot(monocle3_data, aes(x = pseudotime, 
                                                   y = Cluster, 
                                                   color = Cluster)) + 
  geom_point() + theme_classic() + scale_color_manual(name = "Cluster", values = colour_vec) + 
  ggtitle("Pseudotime by cluster", "Monocle 3 Trajectory Analysis")
ggsave("PseudotimeByCluster.pdf", pseudotime_by_cluster, width = 6.5, height = 5.5, units = "in")
ggsave("PseudotimeByCluster.tiff", pseudotime_by_cluster, width = 6.5, height = 5.5, units = "in", dpi = 600)

pseudotime_marker_plot <- plot_cells(monocle_obj, reduction_method = "UMAP",
                                     )


# Slingshot
