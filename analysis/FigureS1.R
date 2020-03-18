library(Seurat)
library(tidyverse)
library(ggpubr)

createSeuratObj <- function(x, sample_name = NULL){
  input_dir <- x
  data <- Read10X(input_dir)
  seurat_obj <- CreateSeuratObject(data, project = "RPE_scRNA")
  
  # Add metadata
  seurat_obj[["condition"]] <- sample_name
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RPS|^RPL")
  return(seurat_obj)
}

calcThreshold <- function(x,
                          nmads = 3,
                          type = c("both", "lower", "upper"),
                          na.rm = FALSE){
  med_val <- stats::median(x, na.rm = na.rm)
  mad_val <- stats::mad(x, center = med_val, na.rm = na.rm)
  upper_limit <- med_val + nmads * mad_val
  lower_limit <- med_val - nmads * mad_val
  
  if (type == "lower"){
    return(lower_limit)
  } else if (type == "upper") {
    return(upper_limit)
  } else(
    return(list(lower_limit, upper_limit))
  )
}

young_obj <- createSeuratObj("/Volumes/LACIE/RPE_scRNA/PUBLISHED/H9_RPE_YOUNG/outs/filtered_feature_bc_matrix/", sample_name = "Young")
aged_obj <- createSeuratObj("/Volumes/LACIE/RPE_scRNA/PUBLISHED/H9_RPE_AGED/outs/filtered_feature_bc_matrix/", sample_name = "Aged")

young_df <- young_obj@meta.data %>% rownames_to_column(var = "cell_barcode")
aged_df <- aged_obj@meta.data %>% rownames_to_column(var = "cell_barcode")

young_df$cell_barcode <- paste0(young_df$cell_barcode, "_1")
aged_df$cell_barcode <- paste0(aged_df$cell_barcode, "_2")

nCount_lower_threshold1 <- calcThreshold(log10(young_df$nCount_RNA), 
                                        nmads = 3, type = "lower")
nCount_upper_threshold1 <- calcThreshold(log10(young_df$nCount_RNA), 
                                        nmads = 3, type = "upper")

nCount_lower_threshold2 <- calcThreshold(log10(aged_df$nCount_RNA), 
                                         nmads = 3, type = "lower")
nCount_upper_threshold2 <- calcThreshold(log10(aged_df$nCount_RNA), 
                                         nmads = 3, type = "upper")

# Convert to count
nCount_lower_threshold1 <- 10^(nCount_lower_threshold1)
nCount_lower_threshold2 <- 10^(nCount_lower_threshold2)
nCount_upper_threshold1 <- 10^(nCount_upper_threshold1)
nCount_upper_threshold2 <- 10^(nCount_upper_threshold2)

umi_plot1 <- ggplot(young_df, aes(nCount_RNA)) + 
  geom_histogram(binwidth = 250, color = "black", fill = "gray23") + 
  theme_classic() + ggtitle("1-month timepoint", subtitle = "Total UMI count per cell") +
  geom_vline(xintercept = nCount_lower_threshold1, color = "#CC0000") + 
  geom_vline(xintercept = nCount_upper_threshold1, color = "#CC0000") +
  xlab("Total UMI counts") + ylab("Number of cells")

umi_plot2 <- ggplot(aged_df, aes(nCount_RNA)) + 
  geom_histogram(binwidth = 250, color = "black", fill = "gray23") + 
  theme_classic() + ggtitle("12-month timepoint", subtitle = "Total UMI count per cell") +
  geom_vline(xintercept = nCount_lower_threshold2, color = "#CC0000") + 
  geom_vline(xintercept = nCount_upper_threshold2, color = "#CC0000") +
  xlab("Total UMI counts") + ylab("Number of cells")

feature_plot1 <- ggplot(young_df, aes(nFeature_RNA)) + 
  geom_histogram(binwidth = 100, color = "black", fill = "gray23") + 
  theme_classic() + ggtitle("1-month timepoint", subtitle = "Number of detected features per cell") +
  geom_vline(xintercept = 220, color = "#CC0000") + 
  geom_vline(xintercept = 5000, color = "#CC0000") +
  xlab("Number of detected features") + ylab("Number of cells")

feature_plot2 <- ggplot(aged_df, aes(nFeature_RNA)) + 
  geom_histogram(binwidth = 100, color = "black", fill = "gray23") + 
  theme_classic() + ggtitle("12-month timepoint", subtitle = "Number of detected features per cell") +
  geom_vline(xintercept = 220, color = "#CC0000") + 
  xlab("Number of detected features") + ylab("Number of cells")

mt_plot1 <- ggplot(young_df, aes(percent.mt)) + 
  geom_histogram(binwidth = 1, color = "black", fill = "gray23") + 
  theme_classic() + ggtitle("1-month timepoint", subtitle = "Proportion of mitochondrial expression to total expression") +
  geom_vline(xintercept = 25, color = "#CC0000") + 
  xlab("% Mt-associated expression") + ylab("Number of cells")

mt_plot2 <- ggplot(aged_df, aes(percent.mt)) + 
  geom_histogram(binwidth = 1, color = "black", fill = "gray23") + 
  theme_classic() + ggtitle("12-month timepoint", subtitle = "Proportion of mitochondrial expression to total expression") +
  geom_vline(xintercept = 25, color = "#CC0000") + 
  xlab("% Mt-associated expression") + ylab("Number of cells")

rb_plot1 <- ggplot(young_df, aes(percent.rb)) + 
  geom_histogram(binwidth = 1, color = "black", fill = "gray23") + 
  theme_classic() + ggtitle("1-month timepoint", subtitle = "Proportion of ribosomal expression to total expression") +
  geom_vline(xintercept = 60, color = "#CC0000") + 
  xlab("% Rb-associated expression") + ylab("Number of cells")

rb_plot2 <- ggplot(aged_df, aes(percent.rb)) + 
  geom_histogram(binwidth = 1, color = "black", fill = "gray23") + 
  theme_classic() + ggtitle("12-month timepoint", subtitle = "Proportion of ribosomal expression to total expression") +
  geom_vline(xintercept = 60, color = "#CC0000") + 
  xlab("% Rb-associated expression") + ylab("Number of cells")

count_feature_plot1 <- ggscatter(young_df, x = "nCount_RNA", 
                                y = "nFeature_RNA", 
                                color = "percent.mt", size = 1, 
                                title = "1-month timepoint",
                                xlab = "Total UMI count",
                                ylab = "Total number of features") + 
  geom_hline(yintercept = 220, color = "#94221f") + 
  geom_hline(yintercept = 5000, color = "#94221f") + 
  geom_vline(xintercept = nCount_lower_threshold1, color = "#94221f") + 
  geom_vline(xintercept = nCount_upper_threshold1, color = "#94221f") +
  gradient_color(c("#1565C0", "#b92b27"))
count_feature_plot1 <- ggpar(count_feature_plot1, 
                            legend = "bottom", legend.title = "% Mt Expression")
count_feature_plot1 <- ggMarginal(count_feature_plot1, type = "histogram")

count_feature_plot2 <- ggscatter(aged_df, x = "nCount_RNA", 
                                 y = "nFeature_RNA", 
                                 color = "percent.mt", size = 1, 
                                 title = "12-month timepoint",
                                 xlab = "Total UMI count",
                                 ylab = "Total number of features") + 
  geom_hline(yintercept = 220, color = "#94221f") + 
  geom_vline(xintercept = nCount_lower_threshold2, color = "#94221f") + 
  geom_vline(xintercept = nCount_upper_threshold2, color = "#94221f") +
  gradient_color(c("#1565C0", "#b92b27"))
count_feature_plot2 <- ggpar(count_feature_plot2, 
                             legend = "bottom", legend.title = "% Mt Expression")
count_feature_plot2 <- ggMarginal(count_feature_plot2, type = "histogram")

tiff("FigureS1.tiff", width = 340, height = 360, units = "mm", res = 600)
ggarrange(ggarrange(umi_plot1, umi_plot2, ncol = 2), 
          ggarrange(feature_plot1, feature_plot2, ncol = 2), 
          ggarrange(mt_plot1, mt_plot2, ncol = 2), 
          ggarrange(rb_plot1, rb_plot2, ncol = 2), 
          ggarrange(count_feature_plot1, count_feature_plot2, ncol = 2), 
          nrow = 5, labels = c("A", "B", "C", "D", "E"), heights = c(1, 1, 1, 1, 2),
          font.label = list(size = 24))
dev.off()
