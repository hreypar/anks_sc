library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(viridis)
library(ComplexHeatmap)
library(ggpubr)


integrated <- readRDS("results/b_integration_azimuth/integrated_seurat_harmony.rds")

condition_colors <- c("Healthy" = "#0173B2", "AS" = "#DE8F05")

theme_publication <- theme_classic() +
  theme(
    text = element_text(size = 12),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 11, face = "bold"),
    strip.background = element_rect(fill = "grey95", color = "grey20"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5)
  )

# SUBSET AND RECLUSTER MONOCYTES
create_monocyte_subclusters <- function(seurat_obj) {
  # Subset ONLY CD14 and CD16 monocytes (exclude cDC2)
  mono_subset <- subset(seurat_obj, 
                        predicted.celltype.l2 %in% c("CD14 Mono", "CD16 Mono"))
  
  # Recalculate variable features focusing on this subset
  mono_subset <- NormalizeData(mono_subset)
  mono_subset <- FindVariableFeatures(mono_subset, selection.method = "vst", nfeatures = 2000)
  
  # ensure S100A8 and S100A9 are included in variable features
  var_features <- unique(c(VariableFeatures(mono_subset), "S100A8", "S100A9"))
  VariableFeatures(mono_subset) <- var_features
  
  mono_subset <- ScaleData(mono_subset)
  mono_subset <- RunPCA(mono_subset, features = VariableFeatures(mono_subset))
  
  mono_subset <- RunUMAP(mono_subset, dims = 1:20)
  mono_subset <- FindNeighbors(mono_subset, dims = 1:20)
  
  # try multiple resolutions to find optimal clustering
  for(res in c(0.3, 0.5, 0.8, 1.0)) {
    mono_subset <- FindClusters(mono_subset, resolution = res)
  }
  
  return(mono_subset)
}

mono_subset <- create_monocyte_subclusters(integrated)

#  compare resolutions
p_res_comparison <- (
  DimPlot(mono_subset, group.by = "RNA_snn_res.0.3", label = TRUE) + 
    ggtitle("Resolution 0.3") + NoLegend() |
    DimPlot(mono_subset, group.by = "RNA_snn_res.0.5", label = TRUE) + 
    ggtitle("Resolution 0.5") + NoLegend() |
    DimPlot(mono_subset, group.by = "RNA_snn_res.0.8", label = TRUE) + 
    ggtitle("Resolution 0.8") + NoLegend() |
    DimPlot(mono_subset, group.by = "RNA_snn_res.1", label = TRUE) + 
    ggtitle("Resolution 1.0") + NoLegend()
)
ggsave("Monocyte_Clustering_Resolution_Comparison.png", p_res_comparison, width = 16, height = 4, dpi = 300)
ggsave("Monocyte_Clustering_Resolution_Comparison.pdf", p_res_comparison, width = 16, height = 4, dpi = 300)


Idents(mono_subset) <- "RNA_snn_res.0.8"

# create a new column to make it clear these are monocyte-specific clusters
mono_subset$monocyte_clusters <- Idents(mono_subset)

table(Idents(mono_subset))


# Discretize S100s EXPRESSION PATTERNS
# Define expression groups
mono_subset$S100A8_group <- cut(
  FetchData(mono_subset, "S100A8")$S100A8,
  breaks = quantile(FetchData(mono_subset, "S100A8")$S100A8, probs = c(0, 0.25, 0.75, 1)),
  labels = c("Low", "Medium", "High"),
  include.lowest = TRUE
)

mono_subset$S100A9_group <- cut(
  FetchData(mono_subset, "S100A9")$S100A9,
  breaks = quantile(FetchData(mono_subset, "S100A9")$S100A9, probs = c(0, 0.25, 0.75, 1)),
  labels = c("Low", "Medium", "High"),
  include.lowest = TRUE
)


# visualize monocytes
p1 <- DimPlot(mono_subset, group.by = "monocyte_clusters", label = TRUE) +
  ggtitle("Monocyte Subclusters") +
  theme_minimal()

p2 <- DimPlot(mono_subset, group.by = "predicted.celltype.l2") +
  ggtitle("Original Cell Type") +
  scale_color_manual(values = c("CD14 Mono" = "#FF6B6B", "CD16 Mono" = "#4ECDC4")) +
  theme_minimal()

p3 <- FeaturePlot(mono_subset, features = "S100A8", order = TRUE) +
  scale_color_viridis(option = "inferno") +
  ggtitle("S100A8 Expression")

p4 <- FeaturePlot(mono_subset, features = "S100A9", order = TRUE) +
  scale_color_viridis(option = "inferno") +
  ggtitle("S100A9 Expression")

p5 <- DimPlot(mono_subset, group.by = "condition", cols = condition_colors) +
  ggtitle("Condition") +
  theme_minimal()

p6 <- DimPlot(mono_subset, group.by = "S100A8_group", 
              cols = c("Low" = "#E6F3FF", "Medium" = "#66A3FF", "High" = "#0040FF")) +
  ggtitle("S100A8 Expression Groups")

# S100A9 expression groups plot
p7 <- DimPlot(mono_subset, group.by = "S100A9_group", 
              cols = c("Low" = "#FFE6E6", "Medium" = "#FF6666", "High" = "#CC0000")) +
  ggtitle("S100A9 Expression Groups")

# Combine plots 
umap_combined <- (p1 | p2 | p5) / (p3 | p4) / (p6 | p7)
ggsave("Monocyte_Subcluster_UMAP_Overview.pdf", umap_combined, width = 15, height = 15, dpi = 300)
ggsave("Monocyte_Subcluster_UMAP_Overview.png", umap_combined, width = 15, height = 15, dpi = 300)


# output them individually 
ggsave("Monocyte_Subclusters.pdf", p1, width = 7, height = 6, dpi = 300)
ggsave("Monocyte_Subclusters.png", p1, width = 7, height = 6, dpi = 300)

ggsave("Monocyte_Original_CellType.pdf", p2, width = 7, height = 6, dpi = 300)
ggsave("Monocyte_Original_CellType.png", p2, width = 7, height = 6, dpi = 300)

ggsave("Monocyte_S100A8_Expression.pdf", p3, width = 7, height = 6, dpi = 300)
ggsave("Monocyte_S100A8_Expression.png", p3, width = 7, height = 6, dpi = 300)

ggsave("Monocyte_S100A9_Expression.pdf", p4, width = 7, height = 6, dpi = 300)
ggsave("Monocyte_S100A9_Expression.png", p4, width = 7, height = 6, dpi = 300)

ggsave("Monocyte_Condition.pdf", p5, width = 7, height = 6, dpi = 300)
ggsave("Monocyte_Condition.png", p5, width = 7, height = 6, dpi = 300)

ggsave("Monocyte_S100A8_Groups.pdf", p6, width = 7, height = 6, dpi = 300)
ggsave("Monocyte_S100A8_Groups.png", p6, width = 7, height = 6, dpi = 300)

ggsave("Monocyte_S100A9_Groups.pdf", p7, width = 7, height = 6, dpi = 300)
ggsave("Monocyte_S100A9_Groups.png", p7, width = 7, height = 6, dpi = 300)


# clusters comp
cluster_celltype_data <- mono_subset@meta.data %>%
  group_by(monocyte_clusters, predicted.celltype.l2) %>%
  summarise(count = n()) %>%
  group_by(monocyte_clusters) %>%
  mutate(percentage = count / sum(count) * 100)

# Barplot of cluster composition by cell type
p_cluster_celltype <- ggplot(cluster_celltype_data, 
                             aes(x = monocyte_clusters, y = percentage, fill = predicted.celltype.l2)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("CD14 Mono" = "#FF6B6B", "CD16 Mono" = "#4ECDC4")) +
  labs(x = "Cluster", y = "Percentage (%)", fill = "Cell Type",
       title = "Cluster Composition by Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

cluster_condition_data <- mono_subset@meta.data %>%
  group_by(monocyte_clusters, condition) %>%
  summarise(count = n()) %>%
  group_by(monocyte_clusters) %>%
  mutate(percentage = count / sum(count) * 100)

# Barplot of condition distribution by cluster
p_cluster_condition <- ggplot(cluster_condition_data, 
                              aes(x = monocyte_clusters, y = percentage, fill = condition)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = condition_colors) +
  labs(x = "Cluster", y = "Percentage (%)", fill = "Condition",
       title = "Condition Distribution by Cluster") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

# Save 
ggsave("Monocyte_Cluster_CellType_Composition.pdf", p_cluster_celltype, width = 8, height = 6, dpi = 300)
ggsave("Monocyte_Cluster_CellType_Composition.png", p_cluster_celltype, width = 8, height = 6, dpi = 300)

ggsave("Monocyte_Cluster_Condition_Distribution.pdf", p_cluster_condition, width = 8, height = 6, dpi = 300)
ggsave("Monocyte_Cluster_Condition_Distribution.png", p_cluster_condition, width = 8, height = 6, dpi = 300)

combined_barplots <- p_cluster_celltype | p_cluster_condition
ggsave("Monocyte_Cluster_Composition_Combined.pdf", combined_barplots, width = 16, height = 6, dpi = 300)
ggsave("Monocyte_Cluster_Composition_Combined.png", combined_barplots, width = 16, height = 6, dpi = 300)


## Violin plots 
# single violin plot for S100A8
v1 <- VlnPlot(mono_subset, 
              features = "S100A8", 
              group.by = "monocyte_clusters",
              pt.size = 0.1) +
  ggtitle("S100A8 Expression by Monocyte Cluster") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "none")

# single violin plot for S100A9
v2 <- VlnPlot(mono_subset, 
              features = "S100A9", 
              group.by = "monocyte_clusters",
              pt.size = 0.1) +
  ggtitle("S100A9 Expression by Monocyte Cluster") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "none")

ggsave("Monocyte_S100A8_Violin.pdf", v1, width = 8, height = 6, dpi = 300)
ggsave("Monocyte_S100A8_Violin.png", v1, width = 8, height = 6, dpi = 300)

ggsave("Monocyte_S100A9_Violin.pdf", v2, width = 8, height = 6, dpi = 300)
ggsave("Monocyte_S100A9_Violin.png", v2, width = 8, height = 6, dpi = 300)

# side by side comparison
plot_data <- FetchData(mono_subset, vars = c("S100A8", "S100A9", "monocyte_clusters", "predicted.celltype.l2"))

plot_data_long <- plot_data %>%
  pivot_longer(cols = c("S100A8", "S100A9"), 
               names_to = "gene", 
               values_to = "expression")

custom_violin <- ggplot(plot_data_long, aes(x = monocyte_clusters, y = expression, fill = monocyte_clusters)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  facet_wrap(~ gene, scales = "free_y", ncol = 2) +
  theme_minimal() +
  labs(x = "Monocyte Cluster", y = "Expression", title = "S100A8 and S100A9 Expression by Monocyte Cluster") +
  theme(legend.position = "none",
        strip.text = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(angle = 0, hjust = 0.5))

ggsave("Monocyte_S100_Custom_Violin.pdf", custom_violin, width = 12, height = 6, dpi = 300)
ggsave("Monocyte_S100_Custom_Violin.png", custom_violin, width = 12, height = 6, dpi = 300)


# S100 correlations by cluster

analyze_s100_correlation_by_cluster <- function(seu) {
  # Get expression data
  expr_data <- FetchData(seu, 
                         vars = c("S100A8", "S100A9", "monocyte_clusters", 
                                  "predicted.celltype.l2", "condition"))
  
  #correlations by cluster and condition
  cor_results <- expr_data %>%
    group_by(monocyte_clusters, condition) %>%
    summarise(
      correlation = cor(S100A8, S100A9, method = "spearman"),
      n_cells = n(),
      mean_S100A8 = mean(S100A8),
      mean_S100A9 = mean(S100A9),
      .groups = "drop"
    )
  
  p_scatter <- ggplot(expr_data, aes(x = S100A8, y = S100A9)) +
    geom_point(aes(color = condition), alpha = 0.5, size = 0.8) +
    geom_smooth(aes(color = condition), method = "lm", se = TRUE) +
    facet_wrap(~monocyte_clusters, scales = "free", ncol = 4) +
    scale_color_manual(values = condition_colors) +
    theme_publication +
    labs(title = "S100A8 vs S100A9 Correlation by Monocyte Subcluster")
  
  return(list(plot = p_scatter, stats = cor_results))
}

cor_analysis <- analyze_s100_correlation_by_cluster(mono_subset)
ggsave("S100_Correlation_by_Subcluster.pdf", cor_analysis$plot, width = 12, height = 10, dpi = 300)
ggsave("S100_Correlation_by_Subcluster.png", cor_analysis$plot, width = 12, height = 10, dpi = 300)

#############################################################
##### violin plots
create_cluster_violin <- function(mono_subset) {
  plot_data <- FetchData(mono_subset, 
                         vars = c("S100A8", "S100A9", "monocyte_clusters", "condition"))
  
  plot_data_long <- plot_data %>%
    pivot_longer(cols = c("S100A8", "S100A9"), 
                 names_to = "Gene", 
                 values_to = "Expression")
  
  stats_data <- plot_data_long %>%
    group_by(Gene, monocyte_clusters, condition) %>%
    summarise(
      mean_expr = mean(Expression),
      pct_pos = mean(Expression > 0) * 100,
      n = n(),
      .groups = "drop"
    )
  
  p1 <- ggplot(plot_data_long, 
               aes(x = monocyte_clusters, y = Expression, fill = condition)) +
    geom_violin(scale = "width", adjust = 1.5, trim = TRUE, alpha = 0.8,
                position = position_dodge(width = 0.9)) +
    geom_boxplot(width = 0.15, position = position_dodge(width = 0.9), 
                 outlier.shape = NA, alpha = 0.6) +
    scale_fill_manual(values = condition_colors, name = "Condition") +
    facet_wrap(~Gene, scales = "free_y", ncol = 2) +
    labs(x = "Monocyte Cluster", y = "Expression Level",
         title = "S100A8/A9 Expression Across Monocyte Subclusters") +
    theme_publication +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
    stat_compare_means(aes(group = condition), 
                       method = "wilcox.test", 
                       label = "p.format",
                       size = 3.5,
                       label.y.npc = 0.95)
  
  return(list(plot = p1, stats = stats_data))
}

create_celltype_violin <- function(mono_subset) {
  # Create data for plotting
  plot_data <- FetchData(mono_subset, 
                         vars = c("S100A8", "S100A9", "predicted.celltype.l2", "condition"))
  
  plot_data_long <- plot_data %>%
    pivot_longer(cols = c("S100A8", "S100A9"), 
                 names_to = "Gene", 
                 values_to = "Expression") %>%
    mutate(predicted.celltype.l2 = factor(predicted.celltype.l2, 
                                          levels = c("CD14 Mono", "CD16 Mono")))
  
  stats_data <- plot_data_long %>%
    group_by(Gene, predicted.celltype.l2, condition) %>%
    summarise(
      mean_expr = mean(Expression),
      pct_pos = mean(Expression > 0) * 100,
      n = n(),
      .groups = "drop"
    )
  
  p2 <- ggplot(plot_data_long, 
               aes(x = predicted.celltype.l2, y = Expression, fill = condition)) +
    geom_violin(scale = "width", adjust = 1.5, trim = TRUE, alpha = 0.8,
                position = position_dodge(width = 0.9)) +
    geom_boxplot(width = 0.15, position = position_dodge(width = 0.9), 
                 outlier.shape = NA, alpha = 0.6) +
    scale_fill_manual(values = condition_colors, name = "Condition") +
    facet_wrap(~Gene, scales = "free_y", ncol = 2) +
    labs(x = "Cell Type", y = "Expression Level",
         title = "S100A8/A9 Expression by Original Cell Type") +
    theme_publication +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    stat_compare_means(aes(group = condition), 
                       method = "wilcox.test", 
                       label = "p.format",
                       size = 3.5,
                       label.y.npc = 0.95)
  
  return(list(plot = p2, stats = stats_data))
}

# Heatmap
create_s100_heatmap <- function(mono_subset) {
  # Calculate mean expression
  heatmap_data <- FetchData(mono_subset, 
                            vars = c("S100A8", "S100A9", "monocyte_clusters", "condition")) %>%
    group_by(monocyte_clusters, condition) %>%
    summarise(
      S100A8_mean = mean(S100A8),
      S100A9_mean = mean(S100A9),
      .groups = "drop"
    ) %>%
    pivot_longer(cols = c(S100A8_mean, S100A9_mean),
                 names_to = "Gene",
                 values_to = "Mean_Expression") %>%
    mutate(Gene = gsub("_mean", "", Gene))
  
  # Create heatmap
  p4 <- ggplot(heatmap_data, 
               aes(x = interaction(condition, monocyte_clusters), 
                   y = Gene, 
                   fill = Mean_Expression)) +
    geom_tile(color = "white", size = 0.5) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = median(heatmap_data$Mean_Expression)) +
    labs(x = "Condition.Cluster", y = "",
         title = "Mean S100A8/A9 Expression",
         fill = "Expression") +
    theme_publication +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p4)
}

# boxplts clusters
create_cluster_comparison <- function(mono_subset) {
  plot_data <- FetchData(mono_subset, 
                         vars = c("S100A8", "S100A9", "monocyte_clusters", 
                                  "condition", "predicted.celltype.l2"))
  
  # combined S100 score
  plot_data$S100_combined <- (plot_data$S100A8 + plot_data$S100A9) / 2
  
  p5 <- ggplot(plot_data, 
               aes(x = monocyte_clusters, y = S100_combined, fill = condition)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(position = position_dodge(width = 0.75), alpha = 0.3, size = 0.5) +
    scale_fill_manual(values = condition_colors) +
    facet_wrap(~predicted.celltype.l2, ncol = 1) +  # Changed to ncol = 1 for vertical stacking
    labs(x = "Monocyte Cluster", y = "Combined S100A8/A9 Expression",
         title = "Combined S100 Expression by Cluster and Cell Type") +
    theme_publication +
    stat_compare_means(aes(group = condition), 
                       method = "wilcox.test", 
                       label = "p.signif",
                       size = 3.5)
  
  return(p5)
}

fig1_cluster_violin <- create_cluster_violin(mono_subset)
fig2_celltype_violin <- create_celltype_violin(mono_subset)
fig4_heatmap <- create_s100_heatmap(mono_subset)
fig5_comparison <- create_cluster_comparison(mono_subset)

ggsave("Monocyte_S100_Cluster_Violin.pdf", fig1_cluster_violin$plot, width = 10, height = 6, dpi = 300)
ggsave("Monocyte_S100_Cluster_Violin.png", fig1_cluster_violin$plot, width = 10, height = 6, dpi = 300)

ggsave("Monocyte_S100_CellType_Violin.pdf", fig2_celltype_violin$plot, width = 8, height = 6, dpi = 300)
ggsave("Monocyte_S100_CellType_Violin.png", fig2_celltype_violin$plot, width = 8, height = 6, dpi = 300)

ggsave("Monocyte_S100_Correlation.pdf", fig3_correlation$plot, width = 12, height = 10, dpi = 300)
ggsave("Monocyte_S100_Correlation.png", fig3_correlation$plot, width = 12, height = 10, dpi = 300)

ggsave("Monocyte_S100_Heatmap.pdf", fig4_heatmap, width = 8, height = 4, dpi = 300)
ggsave("Monocyte_S100_Heatmap.png", fig4_heatmap, width = 8, height = 4, dpi = 300)

ggsave("Monocyte_S100_Combined_Comparison.pdf", fig5_comparison, width = 10, height = 6, dpi = 300)
ggsave("Monocyte_S100_Combined_Comparison.png", fig5_comparison, width = 10, height = 6, dpi = 300)

combined_figure1 <- fig1_cluster_violin$plot / fig3_correlation$plot + 
  plot_layout(heights = c(1, 1.5)) +
  plot_annotation(tag_levels = 'A',
                  title = "S100A8/A9 Expression Analysis in Monocyte Subclusters")

ggsave("Figure_Monocyte_S100_Analysis_Combined.pdf", combined_figure1, width = 12, height = 14, dpi = 300)
ggsave("Figure_Monocyte_S100_Analysis_Combined.png", combined_figure1, width = 12, height = 14, dpi = 300)


create_cluster_cell_counts <- function(mono_subset) {
  # Count cells per cluster
  cell_counts <- mono_subset@meta.data %>%
    group_by(monocyte_clusters) %>%
    summarise(
      n_cells = n(),
      .groups = "drop"
    ) %>%
    mutate(percentage = n_cells / sum(n_cells) * 100)
  
  # Basic barplot
  p_counts <- ggplot(cell_counts, aes(x = monocyte_clusters, y = n_cells)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
    geom_text(aes(label = n_cells), vjust = -0.3, size = 4) +
    labs(x = "Monocyte Cluster", y = "Number of Cells",
         title = "Cell Distribution Across Monocyte Clusters") +
    theme_publication
  
  return(p_counts)
}

# Function to create stacked barplot by condition
create_cluster_condition_counts <- function(mono_subset) {
  # Count cells per cluster and condition
  cell_counts_condition <- mono_subset@meta.data %>%
    group_by(monocyte_clusters, condition) %>%
    summarise(n_cells = n(), .groups = "drop") %>%
    group_by(monocyte_clusters) %>%
    mutate(percentage = n_cells / sum(n_cells) * 100)
  
  # Stacked barplot by condition
  p_stacked <- ggplot(cell_counts_condition, 
                      aes(x = monocyte_clusters, y = n_cells, fill = condition)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = condition_colors) +
    labs(x = "Monocyte Cluster", y = "Number of Cells",
         title = "Cell Distribution by Cluster and Condition") +
    theme_publication
  
  # Add count labels
  label_data <- cell_counts_condition %>%
    group_by(monocyte_clusters) %>%
    arrange(desc(condition)) %>%
    mutate(label_y = cumsum(n_cells) - n_cells/2)
  
  p_stacked <- p_stacked +
    geom_text(data = label_data,
              aes(y = label_y, label = n_cells),
              color = "black", size = 3.5)
  
  return(p_stacked)
}



# Generate all count plots
count_plot1 <- create_cluster_cell_counts(mono_subset)
count_plot2 <- create_cluster_condition_counts(mono_subset)

# Save individual plots
ggsave("Monocyte_Cluster_Cell_Counts.pdf", count_plot1, width = 8, height = 6, dpi = 300)
ggsave("Monocyte_Cluster_Cell_Counts.png", count_plot1, width = 8, height = 6, dpi = 300)

ggsave("Monocyte_Cluster_Condition_Counts.pdf", count_plot2, width = 8, height = 6, dpi = 300)
ggsave("Monocyte_Cluster_Condition_Counts.png", count_plot2, width = 8, height = 6, dpi = 300)


# Combine count plots
combined_counts <- count_plot1 / count_plot2
ggsave("Monocyte_Cluster_Counts_Combined.pdf", combined_counts, width = 10, height = 8, dpi = 300)
ggsave("Monocyte_Cluster_Counts_Combined.png", combined_counts, width = 10, height = 8, dpi = 300)

#################################################################################
mono_subset$condition <- factor(mono_subset$condition, 
                                levels = c("Healthy", "AS"))

# monocytes by condition
p_split <- DimPlot(mono_subset, 
                   split.by = "condition",
                   group.by = "monocyte_clusters",
                   label = TRUE,
                   label.size = 3,
                   repel = TRUE) +
  NoLegend() +
  ggtitle("Monocyte Clusters by Condition")

ggsave("monocyte_umap_split_by_condition.pdf", p_split, width = 14, height = 6, dpi = 300)
ggsave("monocyte_umap_split_by_condition.png", p_split, width = 14, height = 6, dpi = 300)


# highlight specific clusters
clusters_of_interest <- c("0", "5", "6", "7", "11")

for(clust in clusters_of_interest) {
  cells_highlight <- WhichCells(mono_subset, expression = monocyte_clusters == clust)
  
  p_highlight <- DimPlot(mono_subset, 
                         cells.highlight = cells_highlight,
                         split.by = "condition",
                         cols.highlight = "red",
                         cols = "lightgray",
                         pt.size = 0.5) +
    ggtitle(paste("Cluster", clust, "Highlighted"))
  
  ggsave(paste0("monocyte_cluster", clust, "_highlight_by_condition.pdf"), 
         p_highlight, width = 14, height = 6, dpi = 300)
  ggsave(paste0("monocyte_cluster", clust, "_highlight_by_condition.png"), 
         p_highlight, width = 14, height = 6, dpi = 300)
}


# tryensity plots
library(ggpointdensity)

umap_as <- FetchData(subset(mono_subset, condition == "AS"), 
                     vars = c("umap_1", "umap_2"))
umap_healthy <- FetchData(subset(mono_subset, condition == "Healthy"), 
                          vars = c("umap_1", "umap_2"))

p_density_as <- ggplot(umap_as, aes(x = umap_1, y = umap_2)) +
  geom_pointdensity() +
  scale_color_viridis() +
  ggtitle("AS Monocyte Density") +
  theme_minimal() +
  coord_equal() +
  theme(legend.position = "none")

p_density_healthy <- ggplot(umap_healthy, aes(x = umap_1, y = umap_2)) +
  geom_pointdensity() +
  scale_color_viridis() +
  ggtitle("Healthy Monocyte Density") +
  theme_minimal() +
  coord_equal() +
  theme(legend.position = "none")

p_density_combined <- p_density_healthy | p_density_as

ggsave("monocyte_density_by_condition.pdf", p_density_combined, width = 14, height = 6, dpi = 300)
ggsave("monocyte_density_by_condition.png", p_density_combined, width = 14, height = 6, dpi = 300)



#################################################################################
#################################################################################
#
# Markers and DGE
analyze_monocyte_markers <- function(mono_subset) {
  
  # Ensure we are using the monocyte clusters
  Idents(mono_subset) <- "monocyte_clusters"

  all_markers <- FindAllMarkers(mono_subset,
                                only.pos = TRUE,
                                min.pct = 0.25,
                                logfc.threshold = 0.25)
  
  top_markers <- all_markers %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = 20) %>%
    arrange(cluster, desc(avg_log2FC))
  
  # Save marker results
  write.csv(all_markers, "monocyte_cluster_all_markers.csv", row.names = FALSE)
  write.csv(top_markers, "monocyte_cluster_top20_markers.csv", row.names = FALSE)
  
  clusters <- sort(unique(as.character(Idents(mono_subset))))
  pairwise_results <- list()
  
  for(i in 1:(length(clusters)-1)) {
    for(j in (i+1):length(clusters)) {
      clust1 <- clusters[i]
      clust2 <- clusters[j]
      comparison_name <- paste0("cluster", clust1, "_vs_cluster", clust2)
      
      message(paste("Comparing cluster", clust1, "vs cluster", clust2))
      
      markers <- FindMarkers(mono_subset,
                             ident.1 = clust1,
                             ident.2 = clust2,
                             min.pct = 0.1,
                             logfc.threshold = 0.25)
      
      markers$gene <- rownames(markers)
      markers$comparison <- comparison_name
      
      pairwise_results[[comparison_name]] <- markers
      
      write.csv(markers, 
                paste0("pairwise_comparison_", comparison_name, ".csv"),
                row.names = TRUE)
    }
  }
  
  all_pairwise <- do.call(rbind, pairwise_results)
  write.csv(all_pairwise, "monocyte_all_pairwise_comparisons.csv", row.names = FALSE)
  

  # heatmap of top markers
  top10_per_cluster <- all_markers %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = 10) %>%
    pull(gene) %>%
    unique()
  
  # scale the genes we need, so we have all the markers
  genes_to_scale <- unique(c(top10_per_cluster, rownames(mono_subset@assays$RNA$scale.data)))
  mono_subset <- ScaleData(mono_subset, features = genes_to_scale)
  
  p_heatmap <- DoHeatmap(mono_subset, 
                         features = top10_per_cluster,
                         group.by = "monocyte_clusters",
                         group.colors = rainbow(length(unique(mono_subset$monocyte_clusters)))) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    theme(axis.text.y = element_text(size = 8))
  
  ggsave("monocyte_cluster_markers_heatmap.pdf", p_heatmap, width = 16, height = 14, dpi = 300)
  ggsave("monocyte_cluster_markers_heatmap.png", p_heatmap, width = 16, height = 14, dpi = 300)
  
  # DotPlot of top markers
  top5_per_cluster <- all_markers %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = 5) %>%
    pull(gene) %>%
    unique()
  
  p_dot <- DotPlot(mono_subset, 
                   features = top5_per_cluster,
                   group.by = "monocyte_clusters") +
    coord_flip() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave("monocyte_cluster_markers_dotplot.pdf", p_dot, width = 10, height = 12, dpi = 300)
  ggsave("monocyte_cluster_markers_dotplot.png", p_dot, width = 10, height = 12, dpi = 300)
  
  # Feature plots for top markers 
  for(clust in unique(all_markers$cluster)) {
    clust_markers <- all_markers %>%
      filter(cluster == clust) %>%
      slice_max(order_by = avg_log2FC, n = 4) %>%
      pull(gene)
    
    p_features <- FeaturePlot(mono_subset, 
                              features = clust_markers,
                              order = TRUE,
                              ncol = 2) & 
      scale_color_viridis(option = "plasma")
    
    # Add a main title to the combined plot
    p_features <- p_features + 
      plot_annotation(title = paste("Cluster", clust, "- Top Marker Genes"),
                      theme = theme(plot.title = element_text(size = 16, 
                                                              face = "bold", 
                                                              hjust = 0.5)))
    
    ggsave(paste0("monocyte_cluster", clust, "_top_markers_featureplot.pdf"), 
           p_features, width = 10, height = 10, dpi = 300)
    ggsave(paste0("monocyte_cluster", clust, "_top_markers_featureplot.png"), 
           p_features, width = 10, height = 10, dpi = 300)
  }
  
  
  
  # just do the first 3 coparisons 
  for(i in 1:min(3, length(pairwise_results))) {
    comp_name <- names(pairwise_results)[i]
    df <- pairwise_results[[comp_name]]
    
    # Create more informative significance categories
    df$sig_category <- case_when(
      df$p_val_adj < 0.05 & df$avg_log2FC > 0.5 ~ "Upregulated",
      df$p_val_adj < 0.05 & df$avg_log2FC < -0.5 ~ "Downregulated",
      TRUE ~ "Not Significant"
    )
    
    # Create the volcano plot with better colors
    p_volcano <- ggplot(df, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
      # Plot all points
      geom_point(aes(color = sig_category), alpha = 0.6, size = 2) +
      scale_color_manual(values = c(
        "Upregulated" = "red",
        "Downregulated" = "blue",
        "Not Significant" = "gray70"
      ), name = "Significance") +
      # Add threshold lines
      geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", alpha = 0.5) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
      theme_minimal() +
      labs(title = gsub("_", " ", comp_name),
           subtitle = paste("Red = upregulated in", strsplit(comp_name, "_vs_")[[1]][1],
                            "| Blue = upregulated in", strsplit(comp_name, "_vs_")[[1]][2]),
           x = "Average Log2 Fold Change",
           y = "-log10(Adjusted P-value)") +
      theme(legend.position = "bottom",
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5, size = 10))
    
    # Select genes to label - top genes from each direction
    top_upregulated <- df %>%
      filter(sig_category == "Upregulated") %>%
      arrange(desc(avg_log2FC)) %>%
      slice_head(n = 5)  # Reduced from 10 to 5 for each direction
    
    top_downregulated <- df %>%
      filter(sig_category == "Downregulated") %>%
      arrange(avg_log2FC) %>%
      slice_head(n = 5)  # Reduced from 10 to 5 for each direction
    
    # Combine top genes from both directions
    genes_to_label <- bind_rows(top_upregulated, top_downregulated)
    
    # Add labels with better positioning
    if(nrow(genes_to_label) > 0) {
      p_volcano <- p_volcano +
        ggrepel::geom_text_repel(
          data = genes_to_label,
          aes(label = gene),
          size = 3,
          max.overlaps = 15,  # Reduced overlaps
          box.padding = 0.5,
          point.padding = 0.3,
          segment.color = "gray50",
          segment.size = 0.3,
          # Push labels away from the center
          nudge_x = ifelse(genes_to_label$avg_log2FC > 0, 0.5, -0.5),
          nudge_y = 0.5
        )
    }
    
    # Save the plot
    ggsave(paste0("volcano_plot_", comp_name, ".pdf"), p_volcano, width = 10, height = 9, dpi = 300)
    ggsave(paste0("volcano_plot_", comp_name, ".png"), p_volcano, width = 10, height = 9, dpi = 300)
  }

    marker_summary <- all_markers %>%
    group_by(cluster) %>%
    summarise(
      n_markers = n(),
      n_strong_markers = sum(avg_log2FC > 1),
      top_marker = gene[which.max(avg_log2FC)],
      top_marker_fc = max(avg_log2FC)
    )
  
  write.csv(marker_summary, "monocyte_cluster_marker_summary.csv", row.names = FALSE)
  
  return(list(
    all_markers = all_markers,
    top_markers = top_markers,
    pairwise_comparisons = pairwise_results,
    marker_summary = marker_summary
  ))
}


mono_markers <- analyze_monocyte_markers(mono_subset)


################################################
# # try analyzing mitochondrial genes specifically 
# mt_genes <- grep("^MT-", rownames(mono_subset), value = TRUE)
# 
# # Create a module score for mitochondrial stress
# mono_subset <- AddModuleScore(mono_subset,
#                               features = list(mt_genes),
#                               name = "MT_score")
# 
# # Visualize
# FeaturePlot(mono_subset, 
#             features = "MT_score1",
#             split.by = "condition") +
#   scale_color_gradient2(low = "blue", mid = "white", high = "red")
# 

