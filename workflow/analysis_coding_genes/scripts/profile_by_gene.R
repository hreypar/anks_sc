profile_gene_AS_HC <- function(seurat_obj, gene_name, output_dir = NULL) {
  
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  
  if(!gene_name %in% rownames(seurat_obj)) {
    stop(paste(gene_name, "not found in dataset"))
  }
  
  if(is.null(output_dir)) {
    output_dir <- paste0("results/b_bygene_analysis/", gene_name, "_AS_HC_analysis_", format(Sys.Date(), "%Y%m%d"))
  }
  dir.create(output_dir, showWarnings = FALSE)
  dir.create(file.path(output_dir, "plots"), showWarnings = FALSE)
  dir.create(file.path(output_dir, "tables"), showWarnings = FALSE)
  
  seurat_obj$condition <- factor(seurat_obj$condition, levels = c("Healthy", "AS"))
  
  cat("Analyzing", gene_name, "in AS vs HC dataset...\n")
  
  expr_data <- FetchData(seurat_obj, 
                         vars = c(gene_name, "condition", "sample_id",
                                  "predicted.celltype.l1", "predicted.celltype.l2", 
                                  "predicted.celltype.l3"))
  
  overall_stats <- expr_data %>%
    group_by(condition) %>%
    summarise(
      mean_expression = mean(.data[[gene_name]]),
      median_expression = median(.data[[gene_name]]),
      sd_expression = sd(.data[[gene_name]]),
      pct_positive = sum(.data[[gene_name]] > 0) / n() * 100,
      n_positive_cells = sum(.data[[gene_name]] > 0),
      n_total_cells = n(),
      .groups = "drop"
    )
  
  write.csv(overall_stats, file.path(output_dir, "tables", "overall_statistics.csv"), row.names = FALSE)
  
  sample_stats <- expr_data %>%
    group_by(sample_id, condition) %>%
    summarise(
      mean_expression = mean(.data[[gene_name]]),
      median_expression = median(.data[[gene_name]]),
      pct_positive = sum(.data[[gene_name]] > 0) / n() * 100,
      n_cells = n(),
      .groups = "drop"
    )
  
  write.csv(sample_stats, file.path(output_dir, "tables", "sample_statistics.csv"), row.names = FALSE)
  
  # Stats by cell type and condition
  celltype_stats <- expr_data %>%
    group_by(predicted.celltype.l2, condition) %>%
    summarise(
      mean_expression = mean(.data[[gene_name]]),
      median_expression = median(.data[[gene_name]]),
      sd_expression = sd(.data[[gene_name]]),
      pct_positive = sum(.data[[gene_name]] > 0) / n() * 100,
      n_positive_cells = sum(.data[[gene_name]] > 0),
      n_total_cells = n(),
      .groups = "drop"
    ) %>%
    filter(n_total_cells >= 20)
  
  write.csv(celltype_stats, file.path(output_dir, "tables", "celltype_statistics.csv"), row.names = FALSE)
  
  #fold changes
  fc_data <- celltype_stats %>%
    select(predicted.celltype.l2, condition, mean_expression, pct_positive) %>%
    pivot_wider(names_from = condition, values_from = c(mean_expression, pct_positive)) %>%
    mutate(
      fold_change_mean = mean_expression_AS / (mean_expression_Healthy + 0.001),
      log2FC_mean = log2(fold_change_mean),
      diff_pct_positive = pct_positive_AS - pct_positive_Healthy
    ) %>%
    arrange(desc(abs(log2FC_mean)))
  
  write.csv(fc_data, file.path(output_dir, "tables", "celltype_fold_changes.csv"), row.names = FALSE)
  

  # Vln plot by cell type
  p1 <- VlnPlot(seurat_obj, 
                features = gene_name,
                group.by = "predicted.celltype.l2",
                split.by = "condition",
                pt.size = 0) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste(gene_name, "Expression by Cell Type in AS vs HC"))
  
  ggsave(file.path(output_dir, "plots", "violin_by_celltype.pdf"), p1, width = 14, height = 6)
  
  top_celltypes <- celltype_stats %>%
    group_by(predicted.celltype.l2) %>%
    summarise(max_pct = max(pct_positive), .groups = "drop") %>%
    filter(max_pct > 10) %>%
    pull(predicted.celltype.l2)
  
  if(length(top_celltypes) > 0) {
    expr_subset <- expr_data %>%
      filter(predicted.celltype.l2 %in% top_celltypes)
    
    p2 <- ggplot(expr_subset, aes(x = predicted.celltype.l2, y = .data[[gene_name]], fill = condition)) +
      geom_boxplot(outlier.alpha = 0.2) +
      scale_fill_manual(values = c("Healthy" = "#3498db", "AS" = "#e74c3c")) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = paste(gene_name, "Expression in Top Expressing Cell Types"),
           x = "Cell Type", y = "Expression Level")
    
    ggsave(file.path(output_dir, "plots", "boxplot_top_celltypes.pdf"), p2, width = 10, height = 6)
  }
  
  # Feature plot
  p3 <- FeaturePlot(seurat_obj, 
                    features = gene_name,
                    split.by = "condition",
                    order = TRUE,
                    min.cutoff = "q10",
                    max.cutoff = "q95")
  
  ggsave(file.path(output_dir, "plots", "feature_plot_split.pdf"), p3, width = 10, height = 5)
  
  # Dot plot
  p4 <- DotPlot(seurat_obj,
                features = gene_name,
                group.by = "predicted.celltype.l2",
                split.by = "condition",
                cols = c("#3498db", "#e74c3c")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste(gene_name, "Expression: AS vs HC"))
  
  ggsave(file.path(output_dir, "plots", "dotplot_split.pdf"), p4, width = 8, height = 10)
  

  # Wilcoxon test by cell type
  test_results <- list()
  
  for(ct in unique(celltype_stats$predicted.celltype.l2)) {
    ct_data <- expr_data %>% filter(predicted.celltype.l2 == ct)
    
    if(sum(ct_data$condition == "AS") >= 10 & sum(ct_data$condition == "Healthy") >= 10) {
      as_expr <- ct_data %>% filter(condition == "AS") %>% pull(gene_name)
      hc_expr <- ct_data %>% filter(condition == "Healthy") %>% pull(gene_name)
      
      test_result <- wilcox.test(as_expr, hc_expr)
      
      test_results[[ct]] <- data.frame(
        celltype = ct,
        p_value = test_result$p.value,
        n_AS = length(as_expr),
        n_HC = length(hc_expr)
      )
    }
  }
  
  if(length(test_results) > 0) {
    test_df <- bind_rows(test_results)
    test_df$p_adjusted = p.adjust(test_df$p_value, method = "BH")
    test_df <- test_df %>% arrange(p_value)
    
    write.csv(test_df, file.path(output_dir, "tables", "statistical_tests.csv"), row.names = FALSE)
  }
  
  sample_celltype_avg <- expr_data %>%
    group_by(sample_id, condition, predicted.celltype.l2) %>%
    summarise(
      mean_expr = mean(.data[[gene_name]]),
      pct_positive = sum(.data[[gene_name]] > 0) / n() * 100,
      n_cells = n(),
      .groups = "drop"
    ) %>%
    filter(n_cells >= 10)
  
  if(length(top_celltypes) > 0 && length(top_celltypes) <= 6) {
    sample_plot_data <- sample_celltype_avg %>%
      filter(predicted.celltype.l2 %in% top_celltypes[1:min(6, length(top_celltypes))])
    
    p5 <- ggplot(sample_plot_data, aes(x = condition, y = mean_expr, color = condition)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.6) +
      facet_wrap(~ predicted.celltype.l2, scales = "free_y", ncol = 3) +
      scale_color_manual(values = c("Healthy" = "#3498db", "AS" = "#e74c3c")) +
      theme_minimal() +
      labs(title = paste(gene_name, "Expression Variability Across Samples"),
           y = "Mean Expression per Sample")
    
    ggsave(file.path(output_dir, "plots", "sample_variability.pdf"), p5, width = 10, height = 8)
  }
  

  # top expressing cell types for this gene
  top_expressing_celltypes <- celltype_stats %>%
    filter(n_total_cells > 100) %>%
    group_by(predicted.celltype.l2) %>%
    summarise(
      mean_pct = mean(pct_positive),
      total_cells = sum(n_total_cells),
      .groups = "drop"
    ) %>%
    filter(mean_pct > 10) %>%
    arrange(desc(mean_pct)) %>%
    head(3) %>%
    pull(predicted.celltype.l2)
  
  # If no cell types meet criteria, take the top one
  if(length(top_expressing_celltypes) == 0) {
    top_expressing_celltypes <- celltype_stats %>%
      filter(n_total_cells > 100) %>%
      group_by(predicted.celltype.l2) %>%
      summarise(mean_pct = mean(pct_positive), .groups = "drop") %>%
      arrange(desc(mean_pct)) %>%
      head(1) %>%
      pull(predicted.celltype.l2)
  }
  
  for(ct in top_expressing_celltypes) {
    ct_safe <- gsub(" ", "_", ct)
    ct_cells <- subset(seurat_obj, predicted.celltype.l2 == ct)
    
    if(ncol(ct_cells) > 100) {
      cat(paste("  Analyzing", ct, "cells...\n"))
      
      Idents(ct_cells) <- "condition"
      
      tryCatch({
        de_genes <- FindMarkers(ct_cells, 
                                ident.1 = "AS",
                                ident.2 = "Healthy",
                                min.pct = 0.1,
                                logfc.threshold = 0.25)
        
        if(gene_name %in% rownames(de_genes)) {
          gene_de_stats <- de_genes[gene_name, , drop = FALSE]
          write.csv(gene_de_stats, 
                    file.path(output_dir, "tables", 
                              paste0(gene_name, "_DE_stats_", ct_safe, ".csv")))
        }
      }, error = function(e) {
        message(paste("DE analysis failed for", ct, ":", e$message))
      })
      
      tryCatch({
        expr_matrix <- GetAssayData(ct_cells, slot = "data")
        target_expr <- expr_matrix[gene_name, ]
        
        if(sd(target_expr) > 0) {
          cor_vals <- apply(expr_matrix, 1, function(x) {
            if(sd(x) > 0) {
              cor(x, target_expr, method = "spearman", use = "complete.obs")
            } else {
              NA
            }
          })
          
          cor_vals <- cor_vals[!is.na(cor_vals)]
          cor_vals <- sort(cor_vals, decreasing = TRUE)
          
          top_correlated <- data.frame(
            gene = names(cor_vals)[1:min(50, length(cor_vals))],
            correlation = cor_vals[1:min(50, length(cor_vals))]
          )
          
          write.csv(top_correlated, 
                    file.path(output_dir, "tables", 
                              paste0("top_correlated_genes_", ct_safe, ".csv")),
                    row.names = FALSE)
        }
      }, error = function(e) {
        message(paste("Correlation analysis failed for", ct, ":", e$message))
      })
    }
  }
  

  if(length(top_expressing_celltypes) > 0) {
    # use the first top expressing cell type for correlation heatmap
    ct <- top_expressing_celltypes[1]
    ct_cells <- subset(seurat_obj, predicted.celltype.l2 == ct)
    
    if(ncol(ct_cells) > 100) {
      expr_matrix <- as.matrix(GetAssayData(ct_cells, slot = "data"))
      target_expr <- expr_matrix[gene_name, ]
      
      cor_vals <- apply(expr_matrix, 1, function(x) {
        if(sd(x) > 0 & sd(target_expr) > 0) {
          cor(x, target_expr, method = "spearman", use = "complete.obs")
        } else {
          NA
        }
      })
      
      cor_vals <- cor_vals[!is.na(cor_vals)]
      cor_vals <- sort(cor_vals, decreasing = TRUE)
      top_genes <- names(cor_vals)[1:min(20, length(cor_vals))]
      
      cor_matrix <- cor(t(expr_matrix[top_genes, ]), method = "spearman", use = "pairwise.complete.obs")
      
      # Plot
      library(corrplot)
      pdf(file.path(output_dir, "plots", "gene_correlation_heatmap.pdf"), width = 10, height = 10)
      corrplot(cor_matrix, 
               method = "color",
               type = "upper",
               order = "hclust",
               tl.cex = 0.8,
               addCoef.col = "black",
               number.cex = 0.7,
               title = paste("Gene correlations in", ct, "cells"),
               mar = c(0,0,2,0))
      dev.off()
    }
  }
  
  # SUMMARY REPORT

  
  report_file <- file.path(output_dir, paste0(gene_name, "_analysis_report.txt"))
  sink(report_file)
  
  cat("GENE EXPRESSION ANALYSIS REPORT\n")
  cat("===============================\n\n")
  cat("Gene:", gene_name, "\n")
  cat("Dataset: Ankylosing Spondylitis vs Healthy Controls\n")
  cat("Analysis Date:", format(Sys.Date(), "%Y-%m-%d"), "\n\n")
  
  cat("DATASET OVERVIEW:\n")
  cat("Total cells:", ncol(seurat_obj), "\n")
  cat("AS cells:", sum(seurat_obj$condition == "AS"), "\n")
  cat("Healthy cells:", sum(seurat_obj$condition == "Healthy"), "\n")
  cat("Number of samples:", length(unique(seurat_obj$sample_id)), "\n\n")
  
  cat("OVERALL EXPRESSION:\n")
  print(overall_stats)
  cat("\n")
  
  cat("TOP EXPRESSING CELL TYPES:\n")
  top_exp <- celltype_stats %>%
    group_by(predicted.celltype.l2) %>%
    summarise(mean_pct = mean(pct_positive), .groups = "drop") %>%
    arrange(desc(mean_pct)) %>%
    head(10)
  print(as.data.frame(top_exp))
  cat("\n")
  
  cat("CELL TYPES ANALYZED FOR DE/CORRELATION:\n")
  if(length(top_expressing_celltypes) > 0) {
    for(ct in top_expressing_celltypes) {
      cat("-", ct, "\n")
    }
  } else {
    cat("No cell types met criteria for detailed analysis\n")
  }
  cat("\n")
  
  cat("LARGEST FOLD CHANGES (AS vs HC):\n")
  print(head(fc_data[, c("predicted.celltype.l2", "log2FC_mean", "diff_pct_positive")], 10))
  
  if(exists("test_df")) {
    cat("\n\nSIGNIFICANT DIFFERENCES (p < 0.05):\n")
    sig_tests <- test_df %>% filter(p_adjusted < 0.05)
    if(nrow(sig_tests) > 0) {
      print(sig_tests)
    } else {
      cat("No significant differences found after multiple testing correction.\n")
    }
  }
  
  sink()
  
  cat("\nAnalysis complete! Results saved to:", output_dir, "\n")
  
  # Return key results
  return(list(
    overall_stats = overall_stats,
    celltype_stats = celltype_stats,
    fold_changes = fc_data,
    statistical_tests = if(exists("test_df")) test_df else NULL,
    top_expressing_celltypes = top_expressing_celltypes
  ))
}


# do it for S100A8
results_S100A8 <- profile_gene_AS_HC(integrated, "S100A8")

# results_S100A9 <- profile_gene_AS_HC(integrated, "S100A9")
# results_TNF <- profile_gene_AS_HC(integrated, "TNF")
# results_RORC <- profile_gene_AS_HC(integrated, "RORC")

key_genes <- c("IL17A", "IL17F", "IL23A", "IL23R", "ERAP1", "ERAP2", 
               "S100A8", "S100A9", "TNF", "RORC", "STAT3", "IL6", "IFNG")

lapply(key_genes, profile_gene_AS_HC, seurat_obj = integrated)

#profile_gene_AS_HC(integrated, "XIST")
