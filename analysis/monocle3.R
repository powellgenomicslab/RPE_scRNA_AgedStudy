library(monocle3)
library(Seurat)
library(tidyverse)

colour_df <- read_csv("~/Documents/Experiments/RPE_scRNA/Analysis/RevisedAnalysis/cluster_colours.csv")
colour_df$Cluster <- factor(colour_df$Cluster, levels = c(paste0("Young_", as.character(0:5)),
                                                          paste0("Aged_", as.character(0:5)),
                                                          paste0("Common_", as.character(1:6)))
                            )
colour_vec <- as.character(colour_df$Colour)
names(colour_vec) <- as.character(colour_df$Cluster)

# Read in analysed object
seurat_obj <- readRDS("/Volumes/LACIE/RPE_scRNA/PUBLISHED/H9_RPE_IntegratedObject.rds")

metadata <- FetchData(seurat_obj, c("condition", "cluster", "seurat_clusters",
                                     "nCount_RNA", "nFeature_RNA", "percent.mt",
                                     "percent.rb"))

metadata$cluster <- factor(metadata$cluster, levels = c(
  paste0("Common", "_", as.character(1:6)),
  paste0("Young", "_", as.character(0:5)),
  paste0("Aged", "_", as.character(0:5))
))

metadata$Months <- metadata$condition
metadata$Months[metadata$Months == "Young"] <- "1"
metadata$Months[metadata$Months == "Aged"] <- "12"
metadata$Months <- factor(metadata$Months, levels = c("1", "12"))

seurat_obj[["cluster"]] <- metadata$cluster
seurat_obj[["Months"]] <- metadata$Months

Idents(seurat_obj) <- "cluster"
colnames(metadata)[2] <- "cluster_identity"


counts <- seurat_obj[["integrated"]]@scale.data
metadata <- metadata[colnames(counts), ]

gene_annotation <- data.frame(gene_short_name = rownames(counts), 
                              num_cells_expressed = Matrix::rowSums(counts > 0), row.names = rownames(counts))

cds <- new_cell_data_set(counts,
                         cell_metadata = metadata,
                         gene_metadata = gene_annotation)

cds <- preprocess_cds(cds, num_dim = 30, norm_method = "none", scaling = FALSE)
plot_pc_variance_explained(cds)

cds <- reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "PCA")
cds <- cluster_cells(cds)
cds <- learn_graph(cds)

cds_coldata <- colData(cds)
cds_coldata$Months <- cds_coldata$condition
cds_coldata$Months[cds_coldata$Months == "Young"] <- 1
cds_coldata$Months[cds_coldata$Months == "Aged"] <- 12
cds_coldata$Months <- factor(cds_coldata$Months, levels = c("1", "12"))
cds_coldata <- DataFrame(cds_coldata)
colData(cds) <- cds_coldata

age_plot <- plot_cells(cds, color_cells_by = "Months",
                       show_trajectory_graph = FALSE,
                       label_cell_groups = FALSE,
                       label_leaves = FALSE,
                       label_branch_points = FALSE,
                       label_roots = FALSE) +
  scale_color_manual(name = "Months", values = c("1" = "#e71837",
                                                 "12" = "#00674f")) +
  ggtitle("Timepoint", subtitle = "Monocle 3 UMAP")

ggsave("Monocle3_Residual_Timepoint_UMAP.tiff", age_plot, 
       width = 160, height = 145, units = "mm", dpi = 150)

proliferative_genes <- c("MKI67", "TOP2A", "PCLAF", "PTTG1", "RRM2", "TPX2", "PBK", "RAX", "CTNNB1",
                         "WNT2B", "WN3", "WNT5A", "WNT7B", "SOX4", "PBM7")
proliferative_genes <- proliferative_genes[proliferative_genes %in% rownames(counts)]
proliferative_monocle_umap <- plot_cells(cds, genes = proliferative_genes,
                                         show_trajectory_graph= FALSE,
                                         label_branch_points = FALSE, 
                                         label_leaves = FALSE,
                                         label_roots = FALSE,
                                         label_cell_groups=FALSE,
                                         label_groups_by_cluster = FALSE)
ggsave("Monocle3_Residual_ProliferativeMarker_UMAP.tiff", proliferative_monocle_umap, width = 210, height = 150, units = "mm", dpi = 150)

cds <- order_cells(cds)

saveRDS(cds,"/Volumes/LACIE/RPE_scRNA/PUBLISHED/H9_RPE_Residual_Monocle3_Object.rds")

cds <- readRDS("/Volumes/LACIE/RPE_scRNA/PUBLISHED/H9_RPE_Residual_Monocle3_Object.rds")

cluster_plot <- plot_cells(cds, color_cells_by = "cluster_identity",
                           show_trajectory_graph = TRUE,
                           label_cell_groups = FALSE,
                           label_leaves = FALSE,
                           label_branch_points = FALSE,
                           label_roots = FALSE) +
  scale_colour_manual(name = "Cluster", values = colour_vec) +
  ggtitle("Clusters", subtitle = "Monocle 3 Trajectory Analysis")
ggsave("Monocle3_Residual_Clusters_UMAP.tiff", cluster_plot, width = 165, height = 145, units = "mm", dpi = 150)


pseudotime_plot <- plot_cells(cds, color_cells_by = "pseudotime",
                           show_trajectory_graph = TRUE,
                           label_cell_groups = FALSE,
                           label_leaves = FALSE,
                           label_branch_points = FALSE,
                           label_roots = FALSE) +
  ggtitle("Clusters", subtitle = "Monocle 3 Trajectory Analysis")
ggsave("Monocle3_Residual_Pseudotime_UMAP.tiff", pseudotime_plot, width = 160, height = 145, units = "mm", dpi = 150)


proliferative_genes <- c("MKI67", "TOP2A", "TPX2", "PTTG1", "PCLAF", "RRM2")
proliferative_monocle_umap <- plot_cells(cds, genes = proliferative_genes,
                                         show_trajectory_graph= TRUE,
                                         label_branch_points = FALSE, 
                                         label_leaves = FALSE,
                                         label_roots = FALSE,
                                         label_cell_groups=FALSE,
                                         label_groups_by_cluster = FALSE)
ggsave("Monocle3_Residual_ProliferativeMarker_UMAP.tiff", proliferative_monocle_umap, 
       width = 160, height = 100, units = "mm", dpi = 300)


rpe_genes <- c("PAX6", "RAX", "SIX3", "MITF", "PMEL", 
               "RGR", "TYR", "RLBP1", "RBP1", "RPE65")
rpe_monocle_umap <- plot_cells(cds, genes = rpe_genes,
                                         show_trajectory_graph= TRUE,
                                         label_branch_points = FALSE, 
                                         label_leaves = FALSE,
                                         label_roots = FALSE,
                                         label_cell_groups=FALSE,
                                         label_groups_by_cluster = FALSE)
ggsave("Monocle3_Residual_RPE_Marker_UMAP.tiff", rpe_monocle_umap, 
       width = 200, height = 150, units = "mm", dpi = 300)

# Generate Pseudotime Progress plot
col_data <- colData(cds)
pseudotime <- pseudotime(cds)
counts <- seurat_obj[["SCT"]]@data

# Add pseudotime to colData
col_data$pseudotime <-pseudotime  

# Generate pseudotime scatter plot
library(ggpubr)
pseudotime_df <- col_data[, c("Months", "cluster_identity", "pseudotime")]
pseudotime_df$cell_barcode <- rownames(pseudotime_df) 
pseudotime_df <- as.data.frame(pseudotime_df)

pseudotime_plot <- ggplot(pseudotime_df, aes(pseudotime, 
                                             fill = cluster_identity)) +
  geom_density() + theme_classic() + 
  scale_fill_manual(name = "Cluster", values = colour_vec) +
  ggtitle("Pseudotime", "Monocle 3 Trajectory Analysis") + facet_grid(rows = "cluster_identity", scales = "free_y") +
  xlab("Pseudotime") + theme(legend.position = "none", strip.text.y = element_text(angle = 0))
ggsave("Monocle3_Pseudotime_TimelinePlot.tiff", pseudotime_plot, width = 130, height = 240, units = "mm", dpi = 150)

pseudotime_plot <- ggplot(pseudotime_df, aes(x = pseudotime, y = cluster_identity, 
                     colour = cluster_identity)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_manual(name = "Cluster", values = colour_vec) + theme_classic() +
  xlab("Monocle 3 Pseudotime") + ylab("Cluster") +
  ggtitle("Cells ordered by Monocle 3 Pseudotime")
ggsave("Monocle3Pseudotime_BeeswarmPlot.png", pseudotime_plot, width = 140, height = 160, units = "mm", dpi = 150)

library(corrplot)


# Identify genes that are differentially expressed across pseudotime

# Moran's I test
# Statistic tells whether cells at similar positions on a trajectory will have
# similar expression levels for a gene being tested
# Intepretation of Moran's I: +1 means that nearby cells will have perfectly
# similar expression, 0 means no correlation, -1 anti-correlated

trajectory_pr_test_df <- graph_test(cds, neighbor_graph="principal_graph", cores=4)
trjaectory_pr_test_df <- trajectory_pr_test_df %>% rownames_to_column(var = "gene_id")
trajectory_pr_de_genes_df <- trjaectory_pr_test_df %>% filter(q_value < 0.05)
trajectory_pr_de_genes_df <- trajectory_pr_de_genes_df %>% 
  arrange(desc(abs(morans_I))) %>%
  dplyr::select(gene_id, gene_short_name, status, everything())
write_tsv(trajectory_pr_de_genes_df, "Monocle3_Residual_Trajectory_DE_Genes.tsv")

trajectory_pr_de_genes_df <- read_tsv("~/Documents/Experiments/RPE_scRNA/Analysis/RevisedAnalysis/trajectory_analysis/Monocle3_Residual_Trajectory_DE_Genes.tsv")



gene_module_df <- find_gene_modules(cds[trajectory_pr_de_genes_df$gene_id,], 
                                    reduction_method = "UMAP", 
                                    random_seed = 1, cores = 3,
                                    resolution = 10^seq(-6, -1))

gene_module_df <- gene_module_df %>% arrange(supermodule, module)
write_tsv(gene_module_df, "Monocle3_Residual_Trajectory_GeneModules.tsv")

cell_group_df <- tibble::tibble(cell = rownames(colData(cds)),
                                cell_group = colData(cds)$cluster_identity)

agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))

agg_mat <- as.matrix(agg_mat)
agg_mat <- scale(agg_mat)

cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
breaks = seq(0, 1, length=101)

library(ComplexHeatmap)
heatmap <- Heatmap(agg_mat, 
                   cluster_rows = function(m) hclust(dist(agg_mat), method = "ward.D2"),
                   row_names_gp = gpar(fontsize = 8),
                   col = cols)
pdf("Monocle3_Residual_TrajectoryDE_GeneModuleHeatmaps.pdf", width = 7, height = 5)
heatmap
dev.off()

# Test if expression is dependent on time
rpe_genes <- c("PAX6", "SIX3", "RAX", "MITF", "PMEL",
               "RPE65", "RLBP1", "RBP1", "TYR", "TYR1", "RGR")

rpe_subset <- cds[rownames(cds) %in% rpe_genes, ]

gene_expression_boxplots <- plot_genes_violin(rpe_subset, group_cells_by="cluster_identity", ncol=5) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + scale_fill_manual(name = "Cluster", values = colour_vec)
ggsave("PseudotimeGeneExpressionViolinPlots.pdf", gene_expression_boxplots, width = 18 , height = 5)

# Check for batch effect
batch_regressed_gene_fits <- fit_models(rpe_subset, model_formula_str = "~cluster_identity + condition", expression_family = "quasipoisson")
batch_regressed_fit_coefs <- coefficient_table(batch_regressed_gene_fits)
batch_regressed_fit_coefs <- batch_regressed_fit_coefs %>% filter(term != "(Intercept)") %>%
  select(gene_short_name, term, q_value, estimate)
fit_test <- evaluate_fits(batch_regressed_gene_fits)

# Compare models
covariate_test <- compare_models(gene_fits, batch_regressed_gene_fits)

# Need to use count data for these
gene_trajectory_plots <- plot_genes_in_pseudotime(rpe_subset, 
                                                 color_cells_by = "cluster_identity", 
                                                 min_expr = 0.5, 
                                                 vertical_jitter = 0.05, horizontal_jitter = 0.05,
                                                 ncol = 2, nrow = 5) +
  scale_color_manual(name = "Cluster", values = colour_vec)

ggsave("GeneTrajectory_Monocle3_Plots.tiff", gene_trajectory_plots, height = 11, width = 8.5, units = "in", dpi = 150)

