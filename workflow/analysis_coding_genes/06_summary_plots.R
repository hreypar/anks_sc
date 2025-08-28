#!/usr/bin/env Rscript

library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(ggpubr)
library(scales)
library(RColorBrewer)
library(gt)
library(gridExtra)

theme_set(theme_cowplot())

pub_colors <- list(
  condition = c("AS" = "#D62728", "Healthy" = "#1F77B4"),
  celltypes = c(colorRampPalette(brewer.pal(12, "Paired"))(30))
)

integrated <- readRDS("results/b_integration_azimuth/integrated_seurat_harmony.rds")

output_dir <- "results/b_summary_plots"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

metadata <- integrated@meta.data

cat("\n=== Dataset Overview ===\n")
cat("Total cells:", nrow(metadata), "\n")
cat("Number of samples:", length(unique(metadata$sample_id)), "\n")
cat("Conditions:", unique(metadata$condition), "\n")
cat("Cell types:", length(unique(metadata$predicted.celltype.l2)), "\n\n")


p1 <- metadata %>%
  group_by(condition) %>%
  summarise(n_cells = n()) %>%
  ggplot(aes(x = condition, y = n_cells, fill = condition)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = comma(n_cells)), vjust = -0.5, size = 6) +
  scale_fill_manual(values = pub_colors$condition) +
  scale_y_continuous(labels = comma, expand = c(0, 0, 0.1, 0)) +
  labs(title = "Total Cells by Condition",
       x = "Condition",
       y = "Number of Cells") +
  theme(legend.position = "none",
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16))

ggsave(file.path(output_dir, "cells_by_condition.png"), p1, width = 8, height = 6, dpi = 300)
ggsave(file.path(output_dir, "cells_by_condition.pdf"), p1, width = 8, height = 6)

sample_counts <- metadata %>%
  group_by(sample_id, condition) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  arrange(condition, desc(n_cells))

p2 <- ggplot(sample_counts, aes(x = reorder(sample_id, n_cells), y = n_cells, fill = condition)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = comma(n_cells)), hjust = -0.2, size = 3) +
  scale_fill_manual(values = pub_colors$condition) +
  scale_y_continuous(labels = comma, expand = c(0, 0, 0.1, 0)) +
  coord_flip() +
  labs(title = "Cells per Sample",
       x = "Sample ID",
       y = "Number of Cells",
       fill = "Condition") +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 16),
        legend.position = "top")

ggsave(file.path(output_dir, "cells_by_sample.png"), p2, width = 10, height = 12, dpi = 300)
ggsave(file.path(output_dir, "cells_by_sample.pdf"), p2, width = 10, height = 12)


celltype_counts <- metadata %>%
  group_by(predicted.celltype.l2, condition) %>%
  summarise(n_cells = n(), .groups = "drop")

p3 <- ggplot(celltype_counts, aes(x = reorder(predicted.celltype.l2, n_cells), y = n_cells, fill = condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = pub_colors$condition) +
  scale_y_continuous(labels = comma) +
  coord_flip() +
  labs(title = "Cell Type Distribution by Condition",
       x = "Cell Type",
       y = "Number of Cells",
       fill = "Condition") +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.position = "top")

ggsave(file.path(output_dir, "celltype_counts_by_condition.png"), p3, width = 12, height = 10, dpi = 300)
ggsave(file.path(output_dir, "celltype_counts_by_condition.pdf"), p3, width = 12, height = 10)

celltype_pct_condition <- metadata %>%
  group_by(condition, predicted.celltype.l2) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(condition) %>%
  mutate(percentage = n / sum(n) * 100) %>%
  ungroup()

p4 <- ggplot(celltype_pct_condition, aes(x = condition, y = percentage, fill = predicted.celltype.l2)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = pub_colors$celltypes) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), expand = c(0, 0)) +
  labs(title = "Cell Type Composition by Condition",
       x = "Condition",
       y = "Percentage",
       fill = "Cell Type") +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.position = "right",
        legend.text = element_text(size = 10))

ggsave(file.path(output_dir, "celltype_percentage_stacked.png"), p4, width = 12, height = 8, dpi = 300)
ggsave(file.path(output_dir, "celltype_percentage_stacked.pdf"), p4, width = 12, height = 8)

p5 <- ggplot(celltype_pct_condition, aes(x = predicted.celltype.l2, y = percentage, fill = condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = pub_colors$condition) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.position = "top") +
  labs(title = "Cell Type Percentages by Condition",
       x = "Cell Type",
       y = "Percentage within Condition",
       fill = "Condition")

ggsave(file.path(output_dir, "celltype_percentage_dodge.png"), p5, width = 14, height = 8, dpi = 300)
ggsave(file.path(output_dir, "celltype_percentage_dodge.pdf"), p5, width = 14, height = 8)


sample_celltype_pct <- metadata %>%
  group_by(sample_id, condition, predicted.celltype.l2) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(sample_id) %>%
  mutate(percentage = n / sum(n) * 100) %>%
  ungroup()

p6 <- ggplot(sample_celltype_pct, aes(x = sample_id, y = percentage, fill = predicted.celltype.l2)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = pub_colors$celltypes) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), expand = c(0, 0)) +
  facet_wrap(~condition, scales = "free_x", ncol = 2) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        strip.text = element_text(size = 14, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 10)) +
  labs(title = "Cell Type Composition by Individual Sample",
       x = "Sample ID",
       y = "Percentage",
       fill = "Cell Type")

ggsave(file.path(output_dir, "celltype_percentage_by_sample.png"), p6, width = 16, height = 10, dpi = 300)
ggsave(file.path(output_dir, "celltype_percentage_by_sample.pdf"), p6, width = 16, height = 10)

p7 <- ggplot(sample_celltype_pct, aes(x = predicted.celltype.l2, y = percentage, fill = condition)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1), alpha = 0.5) +
  scale_fill_manual(values = pub_colors$condition) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.position = "top") +
  labs(title = "Cell Type Percentage Variability Across Samples",
       x = "Cell Type",
       y = "Percentage within Sample",
       fill = "Condition")

ggsave(file.path(output_dir, "celltype_percentage_boxplot.png"), p7, width = 14, height = 8, dpi = 300)
ggsave(file.path(output_dir, "celltype_percentage_boxplot.pdf"), p7, width = 14, height = 8)

# differences in cell type proportions between conditions
stats_results <- data.frame()

for(ct in unique(metadata$predicted.celltype.l2)) {
  # Get percentages for each sample
  ct_data <- sample_celltype_pct %>%
    filter(predicted.celltype.l2 == ct) %>%
    select(sample_id, condition, percentage)
  
  # Add zeros for samples without this cell type
  all_samples <- unique(metadata$sample_id)
  missing_samples <- setdiff(all_samples, ct_data$sample_id)
  
  if(length(missing_samples) > 0) {
    missing_data <- data.frame(
      sample_id = missing_samples,
      condition = metadata$condition[match(missing_samples, metadata$sample_id)],
      percentage = 0
    )
    ct_data <- rbind(ct_data, missing_data)
  }
  
  as_pct <- ct_data$percentage[ct_data$condition == "AS"]
  hc_pct <- ct_data$percentage[ct_data$condition == "Healthy"]
  
  if(length(as_pct) > 2 && length(hc_pct) > 2) {
    test_result <- wilcox.test(as_pct, hc_pct)
    
    stats_results <- rbind(stats_results, data.frame(
      celltype = ct,
      mean_AS = mean(as_pct),
      mean_Healthy = mean(hc_pct),
      median_AS = median(as_pct),
      median_Healthy = median(hc_pct),
      p_value = test_result$p.value,
      n_AS = length(as_pct),
      n_Healthy = length(hc_pct)
    ))
  }
}

stats_results$p_adjusted = p.adjust(stats_results$p_value, method = "BH")
stats_results$significance <- case_when(
  stats_results$p_adjusted < 0.001 ~ "***",
  stats_results$p_adjusted < 0.01 ~ "**",
  stats_results$p_adjusted < 0.05 ~ "*",
  TRUE ~ "ns"
)

write.csv(stats_results %>% arrange(p_adjusted), 
          file.path(output_dir, "celltype_proportion_statistics.csv"),
          row.names = FALSE)


#comprehensive summary
summary_by_celltype <- metadata %>%
  group_by(condition, predicted.celltype.l2) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  group_by(condition) %>%
  mutate(percentage = round(n_cells / sum(n_cells) * 100, 1)) %>%
  ungroup() %>%
  mutate(label = paste0(comma(n_cells), " (", percentage, "%)")) %>%
  select(condition, predicted.celltype.l2, label) %>%
  pivot_wider(names_from = condition, values_from = label, values_fill = "0 (0%)")

if(exists("stats_results")) {
  summary_by_celltype <- summary_by_celltype %>%
    left_join(stats_results %>% select(celltype, p_adjusted, significance), 
              by = c("predicted.celltype.l2" = "celltype")) %>%
    arrange(p_adjusted, predicted.celltype.l2)
}

gt_table <- summary_by_celltype %>%
  gt() %>%
  tab_header(
    title = "Cell Type Distribution by Condition",
    subtitle = "Number of cells (percentage within condition)"
  ) %>%
  cols_label(
    predicted.celltype.l2 = "Cell Type",
    AS = "AS",
    Healthy = "Healthy",
    p_adjusted = "Adjusted p-value",
    significance = "Significance"
  ) %>%
  fmt_number(
    columns = p_adjusted,
    decimals = 4
  ) %>%
  tab_style(
    style = cell_fill(color = "#FFE4E1"),
    locations = cells_body(
      columns = everything(),
      rows = significance %in% c("*", "**", "***")
    )
  )

gtsave(gt_table, file.path(output_dir, "celltype_summary_table_complete.png"))

write.csv(summary_by_celltype, 
          file.path(output_dir, "celltype_summary_complete.csv"),
          row.names = FALSE)

# get positions for significance stars
y_max <- max(sample_celltype_pct$percentage) * 1.1
sig_data <- stats_results %>%
  filter(p_adjusted < 0.05) %>%
  select(celltype, significance, p_adjusted)

p7_with_stars <- ggplot(sample_celltype_pct, aes(x = predicted.celltype.l2, y = percentage, fill = condition)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1), alpha = 0.5) +
  scale_fill_manual(values = pub_colors$condition) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), 
                     limits = c(0, y_max * 1.15)) +  # Extra space for stars
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.position = "top") +
  labs(title = "Cell Type Percentage Variability Across Samples",
       x = "Cell Type",
       y = "Percentage within Sample",
       fill = "Condition")

if(nrow(sig_data) > 0) {
  for(i in 1:nrow(sig_data)) {
    ct <- sig_data$celltype[i]
    sig <- sig_data$significance[i]
    
    x_pos <- which(levels(factor(sample_celltype_pct$predicted.celltype.l2)) == ct)
    
    p7_with_stars <- p7_with_stars +
      annotate("text", x = x_pos, y = y_max, 
               label = sig, size = 6, fontface = "bold")
  }
  
  p7_with_stars <- p7_with_stars +
    labs(caption = "* p < 0.05, ** p < 0.01, *** p < 0.001 (Wilcoxon test with BH correction)")
}

ggsave(file.path(output_dir, "celltype_percentage_boxplot_with_significance.png"), 
       p7_with_stars, width = 14, height = 10, dpi = 300)
ggsave(file.path(output_dir, "celltype_percentage_boxplot_with_significance.pdf"), 
       p7_with_stars, width = 14, height = 10)

if(nrow(sig_data) > 0) {
  # Filter to only significant cell types
  sig_celltypes <- sig_data$celltype
  
  p_sig_only <- sample_celltype_pct %>%
    filter(predicted.celltype.l2 %in% sig_celltypes) %>%
    ggplot(aes(x = predicted.celltype.l2, y = percentage, fill = condition)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 21) +
    geom_point(position = position_jitterdodge(jitter.width = 0.1), alpha = 0.5) +
    scale_fill_manual(values = pub_colors$condition) +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          axis.text.y = element_text(size = 14),
          axis.title = element_text(size = 16),
          legend.position = "top") +
    labs(title = "Cell Types with Significant Proportion Differences",
         x = "Cell Type",
         y = "Percentage within Sample",
         fill = "Condition")
  
  for(i in 1:nrow(sig_data)) {
    ct <- sig_data$celltype[i]
    sig <- sig_data$significance[i]
    p_val <- sig_data$p_adjusted[i]
    
    x_pos <- which(sig_celltypes == ct)
    y_pos <- sample_celltype_pct %>%
      filter(predicted.celltype.l2 == ct) %>%
      pull(percentage) %>%
      max() * 1.1
    
    p_sig_only <- p_sig_only +
      annotate("text", x = x_pos, y = y_pos, 
               label = paste0(sig, "\np=", round(p_val, 4)), 
               size = 4, fontface = "bold", hjust = 0.5)
  }
  
  ggsave(file.path(output_dir, "significant_celltype_differences.png"), 
         p_sig_only, width = 12, height = 8, dpi = 300)
}

summary_table <- metadata %>%
  group_by(condition) %>%
  summarise(
    n_samples = n_distinct(sample_id),
    total_cells = n(),
    mean_cells_per_sample = round(n() / n_distinct(sample_id)),
    .groups = "drop"
  )

for(ct in head(sort(unique(metadata$predicted.celltype.l2)), 5)) {
  ct_counts <- metadata %>%
    filter(predicted.celltype.l2 == ct) %>%
    group_by(condition) %>%
    summarise(n = n(), .groups = "drop")
  
  summary_table[[ct]] <- ct_counts$n[match(summary_table$condition, ct_counts$condition)]
}

table_theme <- ttheme_default(
  core = list(fg_params = list(cex = 1.2)),
  colhead = list(fg_params = list(cex = 1.3, fontface = "bold"))
)

p_table <- tableGrob(summary_table, theme = table_theme, rows = NULL)

ggsave(file.path(output_dir, "summary_table.png"), p_table, width = 12, height = 4, dpi = 300)

combined_plot <- (p1 | p5) / (p3 | p7_with_stars) +
  plot_annotation(
    title = "Single-Cell RNA-seq Dataset Overview",
    subtitle = paste("Total cells:", comma(nrow(metadata)), 
                     "| Samples:", length(unique(metadata$sample_id)),
                     "| Cell types:", length(unique(metadata$predicted.celltype.l2))),
    theme = theme(
      plot.title = element_text(size = 20, face = "bold"),
      plot.subtitle = element_text(size = 14)
    )
  )

ggsave(file.path(output_dir, "combined_summary_figure.png"), 
       combined_plot, width = 20, height = 16, dpi = 300)
ggsave(file.path(output_dir, "combined_summary_figure.pdf"), 
       combined_plot, width = 20, height = 16)

