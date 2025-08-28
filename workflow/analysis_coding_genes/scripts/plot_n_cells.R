library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(RColorBrewer)
library(scales)

# Create output directory
dir.create("results/b_azimuth/plots", showWarnings = FALSE, recursive = TRUE)

extract_celltype_data <- function(seurat_list) {
  all_metadata <- lapply(names(seurat_list), function(sample_name) {
    if (!is.null(seurat_list[[sample_name]])) {
      meta <- seurat_list[[sample_name]]@meta.data
      meta$sample <- sample_name
      return(meta)
    }
  })
  combined_meta <- do.call(rbind, all_metadata)
  return(combined_meta)
}

combined_metadata <- extract_celltype_data(processed_seurat_list)

celltype_summary <- combined_metadata %>%
  group_by(sample, condition, predicted.celltype.l1) %>%
  summarise(n_cells = n(), .groups = "drop")

celltype_summary_l2 <- combined_metadata %>%
  group_by(sample, condition, predicted.celltype.l2) %>%
  summarise(n_cells = n(), .groups = "drop")

unique_celltypes_l1 <- sort(unique(combined_metadata$predicted.celltype.l1))
unique_celltypes_l2 <- sort(unique(combined_metadata$predicted.celltype.l2))
n_colors_l1 <- length(unique_celltypes_l1)
n_colors_l2 <- length(unique_celltypes_l2)

if (n_colors_l1 <= 12) {
  color_palette_l1 <- brewer.pal(max(3, n_colors_l1), "Set3")
} else {
  color_palette_l1 <- c(
    brewer.pal(12, "Set3"),
    brewer.pal(min(8, n_colors_l1 - 12), "Set2"),
    brewer.pal(min(9, max(3, n_colors_l1 - 20)), "Set1")
  )
}
if (length(color_palette_l1) < n_colors_l1) {
  color_palette_l1 <- scales::hue_pal()(n_colors_l1)
}
names(color_palette_l1) <- unique_celltypes_l1

color_palette_l2 <- scales::hue_pal()(n_colors_l2)
names(color_palette_l2) <- unique_celltypes_l2


# counts by celltype L1
p1 <- ggplot(celltype_summary, aes(x = sample, y = n_cells, fill = predicted.celltype.l1)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ condition, scales = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "bottom") +
  labs(title = "Pre-Integration: Cell Type Distribution (Level 1) by Sample",
       x = "Sample", y = "Number of Cells", fill = "Cell Type L1") +
  scale_fill_manual(values = color_palette_l1)

celltype_percentages <- combined_metadata %>%
  group_by(sample, condition, predicted.celltype.l1) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(total_cells = sum(n_cells),
         percentage = (n_cells / total_cells) * 100) %>%
  ungroup()

p3 <- ggplot(celltype_percentages, 
             aes(x = sample, y = percentage, fill = predicted.celltype.l1)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  facet_wrap(~ condition, scales = "free_x", ncol = 2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "bottom",
        panel.grid.major.x = element_blank()) +
  labs(title = "Pre-Integration: Cell Type Proportions (Level 1) by Sample",
       x = "Sample", y = "Percentage (%)", fill = "Cell Type L1") +
  scale_fill_manual(values = color_palette_l1) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  guides(fill = guide_legend(nrow = 2))

summary_stats_l1 <- combined_metadata %>%
  group_by(condition, predicted.celltype.l1) %>%
  summarise(total_cells = n(),
            n_samples = n_distinct(sample),
            mean_cells_per_sample = total_cells / n_samples,
            .groups = "drop") %>%
  mutate(percentage = (total_cells / sum(total_cells)) * 100)

# Plot 4: Overall comparison L1
p4 <- ggplot(summary_stats_l1, aes(x = predicted.celltype.l1, y = percentage, fill = condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Pre-Integration: Overall Cell Type Distribution - Healthy vs AS",
       x = "Cell Type L1", y = "Percentage (%)", fill = "Condition") +
  scale_fill_manual(values = c("Healthy" = "steelblue", "AS" = "coral"))

celltype_percentages_box <- combined_metadata %>%
  group_by(sample, condition, predicted.celltype.l1) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  group_by(sample, condition) %>%
  mutate(percentage = (n_cells / sum(n_cells)) * 100)

p5 <- ggplot(celltype_percentages_box, aes(x = predicted.celltype.l1, y = percentage, fill = condition)) +
  geom_boxplot(alpha = 0.7) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1), alpha = 0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Pre-Integration: Cell Type Proportion Distribution Across Samples",
       x = "Cell Type L1", y = "Percentage (%)", fill = "Condition") +
  scale_fill_manual(values = c("Healthy" = "steelblue", "AS" = "coral"))


#  L2
p2 <- ggplot(celltype_summary_l2, aes(x = sample, y = n_cells, fill = predicted.celltype.l2)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ condition, scales = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "bottom",
        legend.text = element_text(size = 7)) +
  labs(title = "Pre-Integration: Cell Type Distribution (Level 2) by Sample",
       x = "Sample", y = "Number of Cells", fill = "Cell Type L2") +
  scale_fill_manual(values = color_palette_l2) +
  guides(fill = guide_legend(ncol = 3))

celltype_percentages_box_l2 <- combined_metadata %>%
  group_by(sample, condition, predicted.celltype.l2) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  group_by(sample, condition) %>%
  mutate(percentage = (n_cells / sum(n_cells)) * 100)

p6 <- ggplot(celltype_percentages_box_l2, aes(x = predicted.celltype.l2, y = percentage, fill = condition)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.8) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1), alpha = 0.3, size = 0.8) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8),
        legend.position = "top",
        plot.margin = margin(5, 5, 5, 50)) +
  labs(title = "Pre-Integration: Cell Type Proportion Distribution Across Samples (Level 2)",
       x = "Cell Type L2", y = "Percentage (%)", fill = "Condition") +
  scale_fill_manual(values = c("Healthy" = "steelblue", "AS" = "coral"))

major_celltypes_l2 <- celltype_percentages_box_l2 %>%
  group_by(predicted.celltype.l2) %>%
  summarise(mean_percentage = mean(percentage)) %>%
  arrange(desc(mean_percentage)) %>%
  slice_head(n = 15) %>%
  pull(predicted.celltype.l2)

p6_major <- celltype_percentages_box_l2 %>%
  filter(predicted.celltype.l2 %in% major_celltypes_l2) %>%
  ggplot(aes(x = predicted.celltype.l2, y = percentage, fill = condition)) +
  geom_boxplot(alpha = 0.7) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1), alpha = 0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        legend.position = "top") +
  labs(title = "Pre-Integration: Major Cell Type (L2) Proportion Distribution",
       subtitle = "Top 15 most abundant cell types",
       x = "Cell Type L2", y = "Percentage (%)", fill = "Condition") +
  scale_fill_manual(values = c("Healthy" = "steelblue", "AS" = "coral"))

celltype_percentages_l2 <- combined_metadata %>%
  group_by(sample, condition, predicted.celltype.l2) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(total_cells = sum(n_cells),
         percentage = (n_cells / total_cells) * 100) %>%
  ungroup()

summary_stats_l2 <- combined_metadata %>%
  group_by(condition, predicted.celltype.l2) %>%
  summarise(total_cells = n(),
            n_samples = n_distinct(sample),
            mean_cells_per_sample = total_cells / n_samples,
            .groups = "drop") %>%
  mutate(percentage = (total_cells / sum(total_cells)) * 100) %>%
  arrange(desc(percentage))

summary_table <- combined_metadata %>%
  group_by(condition, sample) %>%
  summarise(total_cells = n(), .groups = "drop") %>%
  group_by(condition) %>%
  summarise(n_samples = n(),
            total_cells = sum(total_cells),
            mean_cells_per_sample = mean(total_cells),
            sd_cells_per_sample = sd(total_cells),
            min_cells = min(total_cells),
            max_cells = max(total_cells))


# Save plot objects
saveRDS(p1, file = "results/b_azimuth/plots/pre_integration_celltype_l1_counts_plot.rds")
saveRDS(p2, file = "results/b_azimuth/plots/pre_integration_celltype_l2_counts_plot.rds")
saveRDS(p3, file = "results/b_azimuth/plots/pre_integration_celltype_l1_percentages_plot.rds")
saveRDS(p4, file = "results/b_azimuth/plots/pre_integration_celltype_comparison_conditions_plot.rds")
saveRDS(p5, file = "results/b_azimuth/plots/pre_integration_celltype_l1_boxplot_distribution_plot.rds")
saveRDS(p6, file = "results/b_azimuth/plots/pre_integration_celltype_l2_boxplot_all_plot.rds")
saveRDS(p6_major, file = "results/b_azimuth/plots/pre_integration_celltype_l2_boxplot_major_plot.rds")

# Save as PDF
ggsave("results/b_azimuth/plots/pre_integration_celltype_l1_counts.pdf", p1, width = 12, height = 8)
ggsave("results/b_azimuth/plots/pre_integration_celltype_l2_counts.pdf", p2, width = 14, height = 8)
ggsave("results/b_azimuth/plots/pre_integration_celltype_l1_percentages.pdf", p3, width = 12, height = 8)
ggsave("results/b_azimuth/plots/pre_integration_celltype_comparison_conditions.pdf", p4, width = 10, height = 6)
ggsave("results/b_azimuth/plots/pre_integration_celltype_l1_boxplot_distribution.pdf", p5, width = 10, height = 6)
ggsave("results/b_azimuth/plots/pre_integration_celltype_l2_boxplot_all.pdf", p6, width = 16, height = 8)
ggsave("results/b_azimuth/plots/pre_integration_celltype_l2_boxplot_major.pdf", p6_major, width = 12, height = 8)

# Save as PNG
ggsave("results/b_azimuth/plots/pre_integration_celltype_l1_counts.png", p1, width = 12, height = 8, dpi = 300)
ggsave("results/b_azimuth/plots/pre_integration_celltype_l2_counts.png", p2, width = 14, height = 8, dpi = 300)
ggsave("results/b_azimuth/plots/pre_integration_celltype_l1_percentages.png", p3, width = 12, height = 8, dpi = 300)
ggsave("results/b_azimuth/plots/pre_integration_celltype_comparison_conditions.png", p4, width = 10, height = 6, dpi = 300)
ggsave("results/b_azimuth/plots/pre_integration_celltype_l1_boxplot_distribution.png", p5, width = 10, height = 6, dpi = 300)
ggsave("results/b_azimuth/plots/pre_integration_celltype_l2_boxplot_all.png", p6, width = 16, height = 8, dpi = 300)
ggsave("results/b_azimuth/plots/pre_integration_celltype_l2_boxplot_major.png", p6_major, width = 12, height = 8, dpi = 300)

# Save data
write.csv(celltype_summary, "results/b_azimuth/plots/pre_integration_celltype_l1_counts_data.csv", row.names = FALSE)
write.csv(celltype_summary_l2, "results/b_azimuth/plots/pre_integration_celltype_l2_counts_data.csv", row.names = FALSE)
write.csv(celltype_percentages, "results/b_azimuth/plots/pre_integration_celltype_l1_percentages_data.csv", row.names = FALSE)
write.csv(celltype_percentages_l2, "results/b_azimuth/plots/pre_integration_celltype_l2_percentages_data.csv", row.names = FALSE)
write.csv(summary_stats_l1, "results/b_azimuth/plots/pre_integration_summary_stats_l1.csv", row.names = FALSE)
write.csv(summary_stats_l2, "results/b_azimuth/plots/pre_integration_summary_stats_l2.csv", row.names = FALSE)
write.csv(summary_table, "results/b_azimuth/plots/pre_integration_summary_table_by_condition.csv", row.names = FALSE)

# Save combined metadata
saveRDS(combined_metadata, "results/b_azimuth/plots/pre_integration_combined_metadata.rds")

# Create combined plots
combined_plot_l1 <- (p1 / p3) | p5
combined_plot_l1 <- combined_plot_l1 + 
  plot_annotation(title = "Pre-Integration Cell Type Level 1 Analysis Overview",
                  theme = theme(plot.title = element_text(size = 16, face = "bold")))

combined_plot_l2 <- (p2 / p6_major)
combined_plot_l2 <- combined_plot_l2 + 
  plot_annotation(title = "Pre-Integration Cell Type Level 2 Analysis Overview",
                  theme = theme(plot.title = element_text(size = 16, face = "bold")))

ggsave("results/b_azimuth/plots/pre_integration_combined_overview_l1.pdf", combined_plot_l1, width = 20, height = 12)
ggsave("results/b_azimuth/plots/pre_integration_combined_overview_l1.png", combined_plot_l1, width = 20, height = 12, dpi = 300)
ggsave("results/b_azimuth/plots/pre_integration_combined_overview_l2.pdf", combined_plot_l2, width = 20, height = 12)
ggsave("results/b_azimuth/plots/pre_integration_combined_overview_l2.png", combined_plot_l2, width = 20, height = 12, dpi = 300)

# Print summary
cat("\nPre-Integration Analysis Complete!\n")
cat("Files saved in: results/b_azimuth/plots/\n")
cat("\nSummary Statistics by Condition:\n")
print(summary_table)
cat("\nTop 10 most abundant L2 cell types:\n")
print(head(summary_stats_l2, 10))
