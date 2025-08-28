# FUNCTION TO CREATE VOLCANO PLOT BETWEEN ANY TWO CLUSTERS

create_cluster_volcano <- function(seurat_obj, 
                                   cluster1, 
                                   cluster2,
                                   group_by = "monocyte_clusters",
                                   fc_threshold = 0.5,
                                   pval_threshold = 0.05,
                                   n_label = 10,
                                   label_s100 = TRUE,
                                   save_plot = TRUE,
                                   output_prefix = "volcano") {
  
  cluster1 <- as.character(cluster1)
  cluster2 <- as.character(cluster2)
  
  Idents(seurat_obj) <- group_by
  
  if(!cluster1 %in% levels(Idents(seurat_obj))) {
    stop(paste("Cluster", cluster1, "not found in", group_by))
  }
  if(!cluster2 %in% levels(Idents(seurat_obj))) {
    stop(paste("Cluster", cluster2, "not found in", group_by))
  }
  
  cat(paste("Finding markers between cluster", cluster1, "and cluster", cluster2, "\n"))
  
  markers <- FindMarkers(seurat_obj,
                         ident.1 = cluster1,
                         ident.2 = cluster2,
                         min.pct = 0.1,
                         logfc.threshold = 0.25)
  
  markers$gene <- rownames(markers)
  
  markers$sig_category <- case_when(
    markers$p_val_adj < pval_threshold & markers$avg_log2FC > fc_threshold ~ "Upregulated",
    markers$p_val_adj < pval_threshold & markers$avg_log2FC < -fc_threshold ~ "Downregulated",
    TRUE ~ "Not Significant"
  )
  
  markers$is_s100 <- grepl("^S100", markers$gene)
  
  n_up <- sum(markers$sig_category == "Upregulated")
  n_down <- sum(markers$sig_category == "Downregulated")
  
  p <- ggplot(markers, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    # Nonsignificant points
    geom_point(data = filter(markers, sig_category == "Not Significant"),
               color = "gray80", alpha = 0.5, size = 1.5) +
    # Significant points
    geom_point(data = filter(markers, sig_category != "Not Significant"),
               aes(color = sig_category), alpha = 0.8, size = 2.5) +
    scale_color_manual(values = c(
      "Upregulated" = "#E41A1C",
      "Downregulated" = "#377EB8"
    ), name = "Significance") +

        geom_vline(xintercept = c(-fc_threshold, fc_threshold), 
               linetype = "dashed", alpha = 0.3) +
    geom_hline(yintercept = -log10(pval_threshold), 
               linetype = "dashed", alpha = 0.3) +
    theme_minimal(base_size = 12) +
    labs(
      title = paste("Cluster", cluster1, "vs Cluster", cluster2),
      subtitle = paste0(n_up, " genes up in cluster ", cluster1, 
                        " | ", n_down, " genes up in cluster ", cluster2),
      x = paste0("Log2 Fold Change (", cluster1, "/", cluster2, ")"),
      y = "-log10(Adjusted P-value)"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )
  
  if(label_s100 && any(markers$is_s100)) {
    p <- p + geom_point(data = filter(markers, is_s100),
                        shape = 21, size = 4, stroke = 1.5, 
                        fill = NA, color = "black")
  }
  
  # upregulated
  top_up <- markers %>%
    filter(sig_category == "Upregulated") %>%
    arrange(desc(avg_log2FC * -log10(p_val_adj))) %>%
    slice_head(n = ceiling(n_label/2))
  
  # downregulated
  top_down <- markers %>%
    filter(sig_category == "Downregulated") %>%
    arrange(avg_log2FC * -log10(p_val_adj)) %>%
    slice_head(n = ceiling(n_label/2))
  
  # label significant S100 genes if requested
  if(label_s100) {
    s100_sig <- markers %>%
      filter(is_s100 & sig_category != "Not Significant")
    
    genes_to_label <- bind_rows(top_up, top_down, s100_sig) %>%
      distinct(gene, .keep_all = TRUE)
  } else {
    genes_to_label <- bind_rows(top_up, top_down)
  }
  
  if(nrow(genes_to_label) > 0) {
    p <- p + ggrepel::geom_text_repel(
      data = genes_to_label,
      aes(label = gene, color = sig_category),
      size = 3,
      max.overlaps = 20,
      box.padding = 0.5,
      point.padding = 0.3,
      segment.color = "gray30",
      segment.size = 0.3,
      show.legend = FALSE,
      min.segment.length = 0.1
    )
  }
  
  if(save_plot) {
    filename <- paste0(output_prefix, "_cluster", cluster1, "_vs_", cluster2)
    ggsave(paste0(filename, ".pdf"), p, width = 10, height = 8, dpi = 300)
    ggsave(paste0(filename, ".png"), p, width = 10, height = 8, dpi = 300)
    cat(paste("Saved plots as", filename, "\n"))
  }
  
  return(list(
    plot = p,
    markers = markers,
    n_upregulated = n_up,
    n_downregulated = n_down
  ))
}


# compare clusters 0 and 5
result_0v5 <- create_cluster_volcano(mono_subset, 
                                     cluster1 = 0, 
                                     cluster2 = 5)

# thresholds
result_0v1 <- create_cluster_volcano(mono_subset, 
                                     cluster1 = "0", 
                                     cluster2 = "1",
                                     fc_threshold = 1,      # 2fold change
                                     pval_threshold = 0.01, 
                                     n_label = 20)          # Label more genes

# multiple comparisons
clusters_to_compare <- list(
  c(3, 4),
  c(3, 6),
  c(7, 8),
  c(8, 16)
)

for(pair in clusters_to_compare) {
  create_cluster_volcano(mono_subset, 
                         cluster1 = pair[1], 
                         cluster2 = pair[2],
                         output_prefix = "monocyte_volcano")
}
