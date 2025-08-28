
# THIS IS NOW T CELLS SCRIPT
# REMOVE THE VIOLIN PLOTS FOR MONOCYTES AND YOU ARE SET


library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(ggpubr)
library(scales)
library(viridis)


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

# Define colorblind-friendly palettes
condition_colors <- c("Healthy" = "#0173B2", "AS" = "#DE8F05")


# Define T cell subsets
t_cell_types <- c("CD4 Naive", "CD4 TCM", "CD4 TEM", "CD4 CTL",
                  "CD8 Naive", "CD8 TCM", "CD8 TEM", 
                  "Treg", "MAIT", "gdT")

create_tcell_heatmap <- function(seurat_obj) {
  # Subset T cells
  t_cells <- subset(seurat_obj, predicted.celltype.l2 %in% t_cell_types)
  
  expr_stats <- t_cells@meta.data %>%
    mutate(
      S100A8_expr = FetchData(t_cells, "S100A8")$S100A8,
      S100A9_expr = FetchData(t_cells, "S100A9")$S100A9
    ) %>%
    group_by(predicted.celltype.l2, condition) %>%
    summarise(
      S100A8_mean = mean(S100A8_expr),
      S100A9_mean = mean(S100A9_expr),
      S100A8_pct = mean(S100A8_expr > 0) * 100,
      S100A9_pct = mean(S100A9_expr > 0) * 100,
      n_cells = n()
    )
  
  create_expr_matrix <- function(gene_col, metric) {
    expr_stats %>%
      select(predicted.celltype.l2, condition, !!sym(gene_col)) %>%
      pivot_wider(names_from = condition, values_from = !!sym(gene_col)) %>%
      column_to_rownames("predicted.celltype.l2") %>%
      as.matrix()
  }
  
  mean_mat_s100a8 <- create_expr_matrix("S100A8_mean", "mean")
  mean_mat_s100a9 <- create_expr_matrix("S100A9_mean", "mean")
  
  pct_mat_s100a8 <- create_expr_matrix("S100A8_pct", "pct")
  pct_mat_s100a9 <- create_expr_matrix("S100A9_pct", "pct")
  
  combined_mean <- cbind(mean_mat_s100a8, mean_mat_s100a9)
  colnames(combined_mean) <- c("S100A8_HC", "S100A8_AS", "S100A9_HC", "S100A9_AS")
  
  combined_pct <- cbind(pct_mat_s100a8, pct_mat_s100a9)
  colnames(combined_pct) <- c("S100A8_HC", "S100A8_AS", "S100A9_HC", "S100A9_AS")
  
  col_fun_expr <- colorRamp2(c(0, 0.5, 1, 2), c("white", "#FEE08B", "#F46D43", "#A50026"))
  col_fun_pct <- colorRamp2(c(0, 10, 25, 50), c("white", "#E0F3DB", "#78C679", "#006837"))
  
  ht1 <- Heatmap(combined_mean,
                 name = "Mean\nExpression",
                 col = col_fun_expr,
                 column_title = "Mean Expression Level",
                 row_names_side = "left",
                 column_names_rot = 45,
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   grid.text(sprintf("%.2f", combined_mean[i, j]), x, y, 
                             gp = gpar(fontsize = 8))
                 })
  
  ht2 <- Heatmap(combined_pct,
                 name = "% Expressing",
                 col = col_fun_pct,
                 column_title = "Percentage of Cells Expressing",
                 show_row_names = FALSE,
                 column_names_rot = 45,
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   grid.text(sprintf("%.1f", combined_pct[i, j]), x, y, 
                             gp = gpar(fontsize = 8))
                 })
  
  ht_list <- ht1 + ht2
  
  return(list(heatmap = ht_list, stats = expr_stats))
}

fig2_result <- create_tcell_heatmap(integrated)

pdf("Figure2_Tcell_S100_Heatmap.pdf", width = 12, height = 8)
draw(fig2_result$heatmap, 
     column_title = "S100A8/A9 Expression in T Cell Subsets", 
     column_title_gp = gpar(fontsize = 16, fontface = "bold"))
dev.off()

png("Figure2_Tcell_S100_Heatmap.png", width = 12, height = 8, units = "in", res = 300)
draw(fig2_result$heatmap, 
     column_title = "S100A8/A9 Expression in T Cell Subsets", 
     column_title_gp = gpar(fontsize = 16, fontface = "bold"))
dev.off()



####### T cell violins
create_tcell_violin <- function(seurat_obj) {
  t_cell_types <- c("CD4 Naive", "CD4 TCM", "CD4 TEM", "CD4 CTL",
                    "CD8 Naive", "CD8 TCM", "CD8 TEM", 
                    "Treg", "MAIT", "gdT")
  
  t_cells <- subset(seurat_obj, predicted.celltype.l2 %in% t_cell_types)
  
  plot_data <- FetchData(t_cells, 
                         vars = c("S100A8", "S100A9", "predicted.celltype.l2", "condition"))
  
  plot_data_long <- plot_data %>%
    pivot_longer(cols = c("S100A8", "S100A9"), 
                 names_to = "Gene", 
                 values_to = "Expression") %>%
    mutate(predicted.celltype.l2 = factor(predicted.celltype.l2, 
                                          levels = t_cell_types))
  
  p_tcell_violin <- ggplot(plot_data_long, 
                           aes(x = predicted.celltype.l2, y = Expression, fill = condition)) +
    geom_violin(scale = "width", adjust = 1.5, trim = TRUE, alpha = 0.8,
                position = position_dodge(width = 0.9)) +
    geom_boxplot(width = 0.15, position = position_dodge(width = 0.9), 
                 outlier.shape = NA, alpha = 0.6) +
    scale_fill_manual(values = condition_colors, name = "Condition") +
    facet_wrap(~Gene, scales = "free_y", ncol = 1) +
    labs(x = "T Cell Subset", y = "Expression Level",
         title = "S100A8/A9 Expression in T Cell Subsets") +
    theme_publication +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    stat_compare_means(aes(group = condition), 
                       method = "wilcox.test", 
                       label = "p.signif",
                       size = 3,
                       hide.ns = TRUE,
                       label.y.npc = 0.95)
  
  return(p_tcell_violin)
}

tcell_violin <- create_tcell_violin(integrated)

ggsave("Figure2B_Tcell_Violin_S100.pdf", tcell_violin, width = 12, height = 10, dpi = 300)
ggsave("Figure2B_Tcell_Violin_S100.png", tcell_violin, width = 12, height = 10, dpi = 300)

#######separate plots
create_effector_tcell_violin <- function(seurat_obj) {
  effector_types <- c("CD4 TEM", "CD4 CTL", "CD8 TEM", "CD8 TCM", "MAIT")
  
  t_cells <- subset(seurat_obj, predicted.celltype.l2 %in% effector_types)
  
  plot_data <- FetchData(t_cells, 
                         vars = c("S100A8", "S100A9", "predicted.celltype.l2", "condition"))
  
  plot_data_long <- plot_data %>%
    pivot_longer(cols = c("S100A8", "S100A9"), 
                 names_to = "Gene", 
                 values_to = "Expression") %>%
    mutate(predicted.celltype.l2 = factor(predicted.celltype.l2, 
                                          levels = effector_types))
  
  p <- ggplot(plot_data_long, 
              aes(x = predicted.celltype.l2, y = Expression, fill = condition)) +
    geom_violin(scale = "width", adjust = 1.5, trim = TRUE, alpha = 0.8,
                position = position_dodge(width = 0.9)) +
    geom_boxplot(width = 0.15, position = position_dodge(width = 0.9), 
                 outlier.shape = NA, alpha = 0.6) +
    scale_fill_manual(values = condition_colors, name = "Condition") +
    facet_wrap(~Gene, scales = "free_y", ncol = 2) +
    labs(x = "T Cell Subset", y = "Expression Level",
         title = "S100A8/A9 Expression in Effector/Memory T Cells") +
    theme_publication +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    stat_compare_means(aes(group = condition), 
                       method = "wilcox.test", 
                       label = "p.format",
                       size = 3.5,
                       label.y.npc = 0.95)
  
  return(p)
}

# focus on CD8 TEM and CD4 CT
create_key_tcell_violin <- function(seurat_obj) {
  key_types <- c("CD8 TEM", "CD4 CTL", "CD4 TEM", "Treg")
  
  t_cells <- subset(seurat_obj, predicted.celltype.l2 %in% key_types)
  
  plot_data <- FetchData(t_cells, 
                         vars = c("S100A8", "S100A9", "predicted.celltype.l2", "condition"))
  
  stats_data <- plot_data %>%
    pivot_longer(cols = c("S100A8", "S100A9"), 
                 names_to = "Gene", 
                 values_to = "Expression") %>%
    group_by(Gene, predicted.celltype.l2, condition) %>%
    summarise(
      pct_pos = mean(Expression > 0) * 100,
      mean_expr = mean(Expression),
      median_expr = median(Expression),
      n = n()
    )
  
  plot_data_long <- plot_data %>%
    pivot_longer(cols = c("S100A8", "S100A9"), 
                 names_to = "Gene", 
                 values_to = "Expression") %>%
    mutate(predicted.celltype.l2 = factor(predicted.celltype.l2, 
                                          levels = key_types))
  
  p <- ggplot(plot_data_long, 
              aes(x = predicted.celltype.l2, y = Expression, fill = condition)) +
    geom_violin(scale = "width", adjust = 1.5, trim = TRUE, alpha = 0.8,
                position = position_dodge(width = 0.9)) +
    geom_boxplot(width = 0.15, position = position_dodge(width = 0.9), 
                 outlier.shape = NA, alpha = 0.6) +
    scale_fill_manual(values = condition_colors, name = "Condition") +
    facet_wrap(~Gene, scales = "free_y", ncol = 2) +
    labs(x = "T Cell Subset", y = "Expression Level",
         title = "S100A8/A9 Expression in Key T Cell Populations",
         subtitle = "Focus on effector and regulatory subsets with highest expression") +
    theme_publication +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    stat_compare_means(aes(group = condition), 
                       method = "wilcox.test", 
                       label = "p.format",
                       size = 3.5,
                       label.y.npc = 0.95)
  
  return(list(plot = p, stats = stats_data))
}

# Generate the plots
effector_violin <- create_effector_tcell_violin(integrated)
key_tcell_result <- create_key_tcell_violin(integrated)

# Save the plots
ggsave("Figure2C_Effector_Tcell_Violin_S100.pdf", effector_violin, width = 10, height = 6, dpi = 300)
ggsave("Figure2C_Effector_Tcell_Violin_S100.png", effector_violin, width = 10, height = 6, dpi = 300)

ggsave("Figure2D_Key_Tcell_Violin_S100.pdf", key_tcell_result$plot, width = 10, height = 6, dpi = 300)
ggsave("Figure2D_Key_Tcell_Violin_S100.png", key_tcell_result$plot, width = 10, height = 6, dpi = 300)

