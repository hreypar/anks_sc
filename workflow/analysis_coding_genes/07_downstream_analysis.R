#!/usr/bin/env Rscript

library(Seurat)
library(tidyverse)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(patchwork)
library(viridis)
library(RColorBrewer)
library(ggpubr)
library(pheatmap)
library(cowplot)
library(ComplexHeatmap)
library(circlize)

theme_set(theme_cowplot())

pub_colors <- list(
  condition = c("AS" = "#D62728", "Healthy" = "#1F77B4"),
  celltypes = c(colorRampPalette(brewer.pal(12, "Paired"))(30)),
  heatmap = colorRampPalette(c("#053061", "#2166AC", "#4393C3", 
                               "#92C5DE", "#D1E5F0", "#FFFFFF", 
                               "#FDDBC7", "#F4A582", "#D6604D", 
                               "#B2182B", "#67001F"))(100)
)

integrated <- readRDS("results/b_integration_azimuth/integrated_seurat_harmony.rds")

if(!"predicted.celltype.l2" %in% colnames(integrated@meta.data)) {
  stop("Cell type annotations not found!")
}

output_dir <- "results/c_differential_expression_comprehensive"
plot_objects_dir <- file.path(output_dir, "plot_objects")
pathway_dir <- file.path(output_dir, "pathway_analysis")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_objects_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pathway_dir, recursive = TRUE, showWarnings = FALSE)

available_reductions <- names(integrated@reductions)
umap_reduction <- if("umap" %in% available_reductions) "umap" else available_reductions[1]

print(table(integrated$condition))

DefaultAssay(integrated) <- "RNA"

if(!"data" %in% Layers(integrated, assay = "RNA")) {
  integrated <- NormalizeData(integrated)
}

integrated$condition <- factor(integrated$condition, 
                               levels = c("Healthy", "AS"))

genes_of_interest <- list(
  IL17_IL23 = c("IL17A", "IL17F", "IL17C", "IL23A", "IL23R", "IL21", "IL22", 
                "RORC", "STAT3", "IL6", "IL6R", "IL1B", "IL12B"),
  Antigen_Processing = c("ERAP1", "ERAP2", "TAP1", "TAP2", "TAPBP", "B2M", "HLA-B"),
  Innate_Immunity = c("S100A8", "S100A9", "S100A12", "TLR4", "TLR2", "NLRP3", "CARD9", "TREM1"),
  TNF_Signaling = c("TNF", "TNFRSF1A", "TNFRSF1B", "NFKB1", "RELA"),
  Interferon = c("IFNA1", "IFNB1", "IRF7", "MX1", "ISG15", "RSAD2"),
  T_Cell = c("TBX21", "GATA3", "FOXP3", "IL10", "TGFB1", "IFNG", "IL4"),
  Bone_Remodeling = c("TNFSF11", "TNFRSF11B", "DKK1", "SOST", "BMP2"),
  Other_AS_Related = c("MEFV", "TNFAIP3", "PTPN22", "CCR1", "CXCR3", "GPR65")
)

key_genes <- c("IL17A", "IL17F", "IL23A", "IL23R", "ERAP1", "ERAP2", 
               "S100A8", "S100A9", "TNF", "RORC", "STAT3", "IL6", "IFNG")

######################### FUNCTIONS

perform_de_analysis <- function(seurat_obj, ident.1 = "AS", ident.2 = "Healthy", 
                                test.use = "wilcox", min.pct = 0.1, 
                                logfc.threshold = 0.25, output_prefix = "all_cells") {
  
  cat(paste("\nPerforming DE analysis for", output_prefix, "...\n"))
  
  n_as <- sum(seurat_obj$condition == ident.1)
  n_hc <- sum(seurat_obj$condition == ident.2)
  cat(paste("Sample sizes -", ident.1, ":", n_as, ",", ident.2, ":", n_hc, "\n"))
  
  Idents(seurat_obj) <- "condition"
  
  # Set logfc.threshold to 0 to get ALL genes
  de_results <- FindMarkers(
    seurat_obj,
    ident.1 = ident.1,
    ident.2 = ident.2,
    test.use = test.use,
    min.pct = min.pct,
    logfc.threshold = 0,  # Changed from 0.25 to 0 to include all genes
    verbose = TRUE
  )
  
  de_results$gene <- rownames(de_results)
  
  de_results$significance <- case_when(
    de_results$p_val_adj < 0.001 ~ "***",
    de_results$p_val_adj < 0.01 ~ "**",
    de_results$p_val_adj < 0.05 ~ "*",
    TRUE ~ "ns"
  )
  
  # Save FULL results (all genes)
  write.csv(de_results, 
            file.path(output_dir, paste0(output_prefix, "_DE_results_FULL.csv")), 
            row.names = FALSE)
  
  # Also save significant results separately for convenience
  sig_results <- de_results %>% 
    filter(p_val_adj < 0.05)
  
  write.csv(sig_results, 
            file.path(output_dir, paste0(output_prefix, "_DE_results_significant.csv")), 
            row.names = FALSE)
  
  sig_genes <- de_results %>% 
    filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.5)
  
  cat(paste("Total genes tested:", nrow(de_results), "\n"))
  cat(paste("Significant genes (FDR < 0.05):", nrow(sig_results), "\n"))
  cat(paste("Significant genes (FDR < 0.05, |log2FC| > 0.5):", nrow(sig_genes), "\n"))
  
  # Note about extreme p-values (informational only)
  extreme_pvals <- sum(de_results$p_val_adj < 1e-50)
  if(extreme_pvals > 0) {
    cat(paste("Note:", extreme_pvals, "genes have p-values < 1e-50\n"))
  }
  
  # Select top 10 upregulated and top 10 downregulated genes
  top_up_genes <- de_results %>%
    filter(p_val_adj < 0.05, avg_log2FC > 0) %>%
    arrange(desc(avg_log2FC), p_val_adj) %>%
    head(10) %>%
    pull(gene)
  
  top_down_genes <- de_results %>%
    filter(p_val_adj < 0.05, avg_log2FC < 0) %>%
    arrange(avg_log2FC, p_val_adj) %>%
    head(10) %>%
    pull(gene)
  
  top_genes <- c(top_up_genes, top_down_genes)
  
  # Standard volcano plot - NO CAPPING
  volcano <- EnhancedVolcano(
    de_results,
    lab = de_results$gene,
    x = 'avg_log2FC',
    y = 'p_val_adj',
    selectLab = top_genes,
    title = paste('AS vs Healthy -', gsub("_", " ", output_prefix)),
    subtitle = "",
    caption = paste0('n = ', nrow(de_results), ' genes'),
    pCutoff = 0.05,
    FCcutoff = 0.5,
    pointSize = 2.0,
    labSize = 3.0,
    colAlpha = 0.8,
    legendPosition = 'right',
    legendLabSize = 10,
    legendIconSize = 3.0,
    drawConnectors = TRUE,
    widthConnectors = 0.5,
    maxoverlapsConnectors = 20,
    gridlines.major = FALSE,
    gridlines.minor = FALSE
  ) + theme(aspect.ratio = 1)
  
  ggsave(file.path(output_dir, paste0(output_prefix, "_volcano_plot.pdf")),
         volcano, width = 10, height = 10, dpi = 300)
  ggsave(file.path(output_dir, paste0(output_prefix, "_volcano_plot.png")),
         volcano, width = 10, height = 10, dpi = 300)
  
  all_as_genes <- unlist(genes_of_interest)
  
  de_results$AS_gene <- de_results$gene %in% all_as_genes
  de_results$color_group <- case_when(
    de_results$AS_gene & de_results$p_val_adj < 0.05 ~ "AS_sig",
    de_results$AS_gene & de_results$p_val_adj >= 0.05 ~ "AS_nonsig",
    !de_results$AS_gene ~ "Other"
  )
  
  # REORDER: Put AS genes at the end so they plot on top
  de_results_plot <- de_results %>%
    arrange(AS_gene)
  
  keyvals <- c("AS_sig" = "red", "AS_nonsig" = "orange", "Other" = "grey")
  
  top_as_genes <- de_results_plot %>%
    filter(AS_gene, p_val_adj < 0.05) %>%
    arrange(p_val_adj) %>%
    pull(gene)
  
  volcano_as <- EnhancedVolcano(
    de_results_plot,
    lab = de_results_plot$gene,
    x = 'avg_log2FC',
    y = 'p_val_adj',
    selectLab = top_as_genes,
    title = paste('AS vs Healthy -', gsub("_", " ", output_prefix), '\nAS Genes Highlighted'),
    subtitle = "",
    caption = paste0('Red = AS genes (p<0.05), Orange = AS genes (ns)'),
    pCutoff = 0.05,
    FCcutoff = 0.5,
    pointSize = 2.0,
    labSize = 3.5,
    colAlpha = 0.8,
    legendPosition = 'right',
    legendLabSize = 10,
    legendIconSize = 3.0,
    drawConnectors = TRUE,
    boxedLabels = TRUE,
    widthConnectors = 0.5,
    colCustom = keyvals[de_results_plot$color_group],
    gridlines.major = FALSE,
    gridlines.minor = FALSE
  ) + theme(aspect.ratio = 1)
  
  ggsave(file.path(output_dir, paste0(output_prefix, "_volcano_AS_genes_highlighted.pdf")),
         volcano_as, width = 10, height = 10, dpi = 300)
  ggsave(file.path(output_dir, paste0(output_prefix, "_volcano_AS_genes_highlighted.png")),
         volcano_as, width = 10, height = 10, dpi = 300)
  
  # Check AS gene enrichment
  de_as_genes <- de_results %>%
    filter(gene %in% all_as_genes, p_val_adj < 0.05)
  
  cat(paste("\nAS panel genes found significant:", nrow(de_as_genes), 
            "out of", sum(de_results$gene %in% all_as_genes), "tested\n"))
  
  if(nrow(de_as_genes) > 0) {
    write.csv(de_as_genes %>% 
                arrange(p_val_adj),
              file.path(output_dir, paste0(output_prefix, "_AS_genes_DE.csv")),
              row.names = FALSE)
  }
  
  return(de_results)
}

perform_gsea <- function(de_results, output_prefix = "all_cells") {
  
  enrich_dir <- file.path(pathway_dir, output_prefix)
  dir.create(enrich_dir, recursive = TRUE, showWarnings = FALSE)
  
  # GENE SET ENRICHMENT ANALYSIS (GSEA) using ALL genes
  de_results_clean <- de_results %>%
    filter(!is.na(p_val_adj), !is.na(avg_log2FC))
  
  # Create ranking metric: sign(log2FC) * -log10(p-value)
  # This combines both fold change and significance
  de_results_clean$rank_metric <- sign(de_results_clean$avg_log2FC) * -log10(de_results_clean$p_val + 1e-300)
  
  # Alternative ranking: by log2FC only
  de_results_clean$rank_by_fc <- de_results_clean$avg_log2FC
  
  # Convert to Entrez IDs
  gene_list_df <- bitr(de_results_clean$gene, 
                       fromType = "SYMBOL", 
                       toType = "ENTREZID", 
                       OrgDb = org.Hs.eg.db)
  
  de_results_entrez <- de_results_clean %>%
    inner_join(gene_list_df, by = c("gene" = "SYMBOL"))
  
  gene_list <- de_results_entrez$rank_metric
  names(gene_list) <- de_results_entrez$ENTREZID
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  gene_list_fc <- de_results_entrez$rank_by_fc
  names(gene_list_fc) <- de_results_entrez$ENTREZID
  gene_list_fc <- sort(gene_list_fc, decreasing = TRUE)
  
  cat(paste("Genes included in GSEA:", length(gene_list), "\n"))
  
  tryCatch({
    gsea_go <- gseGO(geneList = gene_list,
                     OrgDb = org.Hs.eg.db,
                     ont = "BP",
                     minGSSize = 10,
                     maxGSSize = 500,
                     pvalueCutoff = 0.05,
                     verbose = TRUE,
                     eps = 1e-10)
    
    if(!is.null(gsea_go) && nrow(gsea_go@result) > 0) {
      write.csv(gsea_go@result, 
                file.path(enrich_dir, paste0(output_prefix, "_GO_GSEA_results.csv")),
                row.names = FALSE)
      
      # Plot top pathways
      p_gsea <- dotplot(gsea_go, showCategory = 20, split = ".sign") + 
        facet_grid(.~.sign) +
        ggtitle("GO Gene Set Enrichment Analysis")
      ggsave(file.path(enrich_dir, paste0(output_prefix, "_GO_GSEA_dotplot.pdf")), 
             p_gsea, width = 14, height = 10)
      
      # Ridge plot
      p_ridge <- ridgeplot(gsea_go, showCategory = 20) + 
        ggtitle("GO GSEA - Ridge Plot")
      ggsave(file.path(enrich_dir, paste0(output_prefix, "_GO_GSEA_ridgeplot.pdf")), 
             p_ridge, width = 12, height = 10)
      
      # GSEA plots for top pathways
      top_pathways <- gsea_go@result %>%
        arrange(pvalue) %>%
        head(5) %>%
        pull(ID)
      
      for(i in seq_along(top_pathways)) {
        p_gsea_plot <- gseaplot2(gsea_go, geneSetID = top_pathways[i], 
                                 title = gsea_go@result$Description[gsea_go@result$ID == top_pathways[i]])
        ggsave(file.path(enrich_dir, paste0(output_prefix, "_GSEA_plot_", i, ".pdf")), 
               p_gsea_plot, width = 10, height = 8)
      }
    }
  }, error = function(e) {
    cat("Error in GO GSEA:", e$message, "\n")
  })
  
  # KEGG GSEA
  tryCatch({
    gsea_kegg <- gseKEGG(geneList = gene_list,
                         organism = 'hsa',
                         minGSSize = 10,
                         maxGSSize = 500,
                         pvalueCutoff = 0.05,
                         verbose = TRUE)
    
    if(!is.null(gsea_kegg) && nrow(gsea_kegg@result) > 0) {
      write.csv(gsea_kegg@result, 
                file.path(enrich_dir, paste0(output_prefix, "_KEGG_GSEA_results.csv")),
                row.names = FALSE)
      
      p_kegg <- dotplot(gsea_kegg, showCategory = 20, split = ".sign") + 
        facet_grid(.~.sign) +
        ggtitle("KEGG Gene Set Enrichment Analysis")
      ggsave(file.path(enrich_dir, paste0(output_prefix, "_KEGG_GSEA_dotplot.pdf")), 
             p_kegg, width = 14, height = 10)
    }
  }, error = function(e) {
    cat("Error in KEGG GSEA:", e$message, "\n")
  })
  
  # MSigDB GSEA (Hallmark gene sets)
  tryCatch({
    library(msigdbr)
    
    # Get hallmark gene sets
    h_gene_sets <- msigdbr(species = "Homo sapiens", category = "H") %>%
      dplyr::select(gs_name, entrez_gene)
    
    gsea_hallmark <- GSEA(gene_list, TERM2GENE = h_gene_sets, 
                          minGSSize = 10, maxGSSize = 500, 
                          pvalueCutoff = 0.05)
    
    if(!is.null(gsea_hallmark) && nrow(gsea_hallmark@result) > 0) {
      write.csv(gsea_hallmark@result, 
                file.path(enrich_dir, paste0(output_prefix, "_Hallmark_GSEA_results.csv")),
                row.names = FALSE)
      
      p_hallmark <- dotplot(gsea_hallmark, showCategory = 20, split = ".sign") + 
        facet_grid(.~.sign) +
        ggtitle("Hallmark Gene Set Enrichment Analysis")
      ggsave(file.path(enrich_dir, paste0(output_prefix, "_Hallmark_GSEA_dotplot.pdf")), 
             p_hallmark, width = 14, height = 10)
    }
  }, error = function(e) {
    cat("MSigDB analysis skipped (install msigdbr package if needed):", e$message, "\n")
  })
  
  # Save the ranked gene list for external analysis
  ranked_genes_output <- de_results_entrez %>%
    dplyr::select(gene, ENTREZID, avg_log2FC, p_val, p_val_adj, rank_metric) %>%
    arrange(desc(rank_metric))
  
  write.csv(ranked_genes_output,
            file.path(enrich_dir, paste0(output_prefix, "_ranked_genes_for_GSEA.csv")),
            row.names = FALSE)
  
  # Also save in GSEA .rnk format
  rnk_output <- data.frame(
    gene = de_results_clean$gene,
    rank = de_results_clean$rank_metric
  ) %>%
    arrange(desc(rank))
  
  write.table(rnk_output,
              file.path(enrich_dir, paste0(output_prefix, "_GSEA.rnk")),
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  cat(paste("\nGSEA complete. Results saved to:", enrich_dir, "\n"))
  
  return(length(gene_list))
}

analyze_as_pathways <- function(seurat_obj, de_results, output_prefix = "all_cells") {
  
  # module scores for each pathway
  cat("Calculating pathway module scores...\n")
  for(pathway_name in names(genes_of_interest)) {
    pathway_genes <- genes_of_interest[[pathway_name]]
    pathway_genes <- pathway_genes[pathway_genes %in% rownames(seurat_obj)]
    
    if(length(pathway_genes) >= 3) {
      seurat_obj <- AddModuleScore(seurat_obj,
                                   features = list(pathway_genes),
                                   name = paste0(pathway_name, "_Score"))
    }
  }
  
  # AS genes expression data
  all_as_genes <- unlist(genes_of_interest)
  as_genes_present <- all_as_genes[all_as_genes %in% rownames(seurat_obj)]
  
  # get DE results for ALL AS genes
  as_de_summary <- de_results %>%
    filter(gene %in% as_genes_present) %>%
    mutate(
      pathway = case_when(
        gene %in% genes_of_interest$IL17_IL23 ~ "IL17/IL23",
        gene %in% genes_of_interest$Antigen_Processing ~ "Antigen Processing",
        gene %in% genes_of_interest$Innate_Immunity ~ "Innate Immunity",
        gene %in% genes_of_interest$TNF_Signaling ~ "TNF Signaling",
        gene %in% genes_of_interest$Interferon ~ "Interferon",
        gene %in% genes_of_interest$T_Cell ~ "T Cell",
        gene %in% genes_of_interest$Bone_Remodeling ~ "Bone Remodeling",
        gene %in% genes_of_interest$Other_AS_Related ~ "Other AS Related"
      ),
      is_significant = p_val_adj < 0.05
    ) %>%
    arrange(pathway, avg_log2FC)
  
  write.csv(as_de_summary, 
            file.path(pathway_dir, paste0("AS_genes_DE_summary_all_with_pathways_", output_prefix, ".csv")),
            row.names = FALSE)
  
  # create comprehensive barplot with ALL AS genes
  if(nrow(as_de_summary) > 0) {
    pathway_colors <- c(
      "IL17/IL23" = "#E41A1C",
      "Antigen Processing" = "#377EB8", 
      "Innate Immunity" = "#4DAF4A",
      "TNF Signaling" = "#984EA3",
      "Interferon" = "#FF7F00",
      "T Cell" = "#FFFF33",
      "Bone Remodeling" = "#A65628",
      "Other AS Related" = "#F781BF"
    )
    
    p_as_all <- as_de_summary %>%
      ggplot(aes(x = reorder(gene, avg_log2FC), y = avg_log2FC, fill = pathway)) +
      geom_bar(stat = "identity") +
      geom_point(aes(shape = is_significant), size = 3, position = position_nudge(y = 0.02)) +
      scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1), 
                         labels = c("TRUE" = "FDR < 0.05", "FALSE" = "ns")) +
      coord_flip() +
      scale_fill_manual(values = pathway_colors) +
      theme_cowplot() +
      labs(title = paste("All AS Panel Genes - Expression Changes in", gsub("_", " ", output_prefix)),
           subtitle = paste("Total:", nrow(as_de_summary), "genes |", 
                            sum(as_de_summary$is_significant), "significant"),
           x = "Gene",
           y = "Log2 Fold Change (AS vs Healthy)",
           fill = "Pathway",
           shape = "Significance") +
      geom_hline(yintercept = 0, linetype = "dashed") +
      theme(legend.position = "right")
    
    ggsave(file.path(pathway_dir, paste0("AS_genes_all_barplot_", output_prefix, ".pdf")),
           p_as_all, width = 12, height = max(8, nrow(as_de_summary) * 0.25))
  }
  
  # create plot for key genes
  key_genes_de <- as_de_summary %>%
    filter(gene %in% key_genes)
  
  if(nrow(key_genes_de) > 0) {
    p_key_genes <- key_genes_de %>%
      ggplot(aes(x = reorder(gene, avg_log2FC), y = avg_log2FC, fill = pathway)) +
      geom_bar(stat = "identity") +
      geom_point(aes(shape = is_significant), size = 4, position = position_nudge(y = 0.02)) +
      scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1), 
                         labels = c("TRUE" = "FDR < 0.05", "FALSE" = "ns")) +
      coord_flip() +
      scale_fill_manual(values = pathway_colors) +
      theme_cowplot() +
      labs(title = paste("Key AS Genes - Expression Changes in", gsub("_", " ", output_prefix)),
           x = "Gene",
           y = "Log2 Fold Change (AS vs Healthy)",
           fill = "Pathway",
           shape = "Significance") +
      geom_hline(yintercept = 0, linetype = "dashed") +
      theme(legend.position = "right",
            axis.text.y = element_text(size = 12, face = "bold"))
    
    ggsave(file.path(pathway_dir, paste0("key_genes_barplot_", output_prefix, ".pdf")),
           p_key_genes, width = 10, height = 8)
  }
  
  # pathway summary
  cat("\nAS Gene Panel Summary:\n")
  cat("Total AS genes in dataset:", length(as_genes_present), "out of", length(all_as_genes), "\n")
  cat("AS genes differentially expressed:", sum(as_de_summary$p_val_adj < 0.05), "\n")
  cat("\nBreakdown by pathway:\n")
  
  pathway_summary <- data.frame()
  for(pathway_name in names(genes_of_interest)) {
    pathway_genes <- genes_of_interest[[pathway_name]]
    genes_in_data <- sum(pathway_genes %in% rownames(seurat_obj))
    
    pathway_display_name <- case_when(
      pathway_name == "IL17_IL23" ~ "IL17/IL23",
      pathway_name == "Antigen_Processing" ~ "Antigen Processing",
      pathway_name == "Innate_Immunity" ~ "Innate Immunity",
      pathway_name == "TNF_Signaling" ~ "TNF Signaling",
      pathway_name == "Interferon" ~ "Interferon",
      pathway_name == "T_Cell" ~ "T Cell",
      pathway_name == "Bone_Remodeling" ~ "Bone Remodeling",
      pathway_name == "Other_AS_Related" ~ "Other AS Related"
    )
    
    genes_de <- sum(as_de_summary$pathway == pathway_display_name & 
                      as_de_summary$is_significant, na.rm = TRUE)
    
    pathway_summary <- rbind(pathway_summary, data.frame(
      pathway = pathway_name,
      pathway_display = pathway_display_name,
      total_in_panel = length(pathway_genes),
      found_in_data = genes_in_data,
      significant = genes_de,
      percent_sig = ifelse(genes_in_data > 0, round(genes_de/genes_in_data * 100, 1), 0)
    ))
  }
  
  print(pathway_summary %>% dplyr::select(-pathway))
  
  write.csv(pathway_summary, 
            file.path(pathway_dir, paste0("pathway_summary_comprehensive_", output_prefix, ".csv")),
            row.names = FALSE)
  
  pathway_gene_details <- as_de_summary %>%
    dplyr::select(gene, pathway, avg_log2FC, p_val_adj, is_significant) %>%
    arrange(pathway, p_val_adj)
  
  write.csv(pathway_gene_details,
            file.path(pathway_dir, paste0("pathway_genes_detailed_results_", output_prefix, ".csv")),
            row.names = FALSE)
  
  return(seurat_obj)
}

create_individual_feature_plots <- function(seurat_obj, genes, split_by = "condition", 
                                            reduction = umap_reduction) {
  
  genes_present <- genes[genes %in% rownames(seurat_obj)]
  
  if(length(genes_present) == 0) {
    cat("Error: No genes found in the dataset!\n")
    return(NULL)
  }
  
  # create cell type reference split by condition
  celltype_ref <- DimPlot(seurat_obj, 
                          reduction = reduction,
                          group.by = "predicted.celltype.l2",
                          split.by = "condition",
                          label = TRUE,
                          label.size = 3,
                          repel = TRUE,
                          cols = pub_colors$celltypes) +
    ggtitle("Cell Type Distribution: AS vs Healthy") +
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
          legend.position = "none")
  
  ggsave(file.path(output_dir, "celltype_reference_AS_vs_Healthy.pdf"),
         celltype_ref, width = 14, height = 8)
  
  pdf(file.path(output_dir, "feature_plots_key_genes_with_labels.pdf"), width = 14, height = 8)
  
  print(celltype_ref)
  
  # For each gene, create feature plot with cell type labels overlaid
  for(i in seq_along(genes_present)) {
    gene <- genes_present[i]
    cat(paste("Plotting", gene, "(", i, "/", length(genes_present), ")...\n"))
    
    expr_data <- FetchData(seurat_obj, vars = c(gene, "predicted.celltype.l2", "condition", 
                                                paste0(reduction, "_1"), paste0(reduction, "_2")))
    
    # ename gene column to avoid special character issues
    colnames(expr_data)[1] <- "gene_expression"
    
    # cells with expression are plotted on top
    expr_data <- expr_data %>%
      arrange(gene_expression)
    
    p <- ggplot(expr_data, aes_string(x = paste0(reduction, "_1"), 
                                      y = paste0(reduction, "_2"))) +
      geom_point(aes(color = gene_expression), size = 0.5) +
      scale_color_gradient(low = "lightgrey", high = "red", name = "Expression") +
      facet_wrap(~condition) +
      theme_minimal() +
      theme(
        strip.text = element_text(size = 18, face = "bold"),
        legend.position = "right",
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size = 24, face = "bold", hjust = 0.5)
      ) +
      ggtitle(gene)
    
    label_data <- expr_data %>%
      group_by(predicted.celltype.l2, condition) %>%
      summarise(
        x = median(get(paste0(reduction, "_1"))),
        y = median(get(paste0(reduction, "_2"))),
        .groups = "drop"
      )
    
    p <- p + 
      geom_text(data = label_data, 
                aes(x = x, y = y, label = predicted.celltype.l2),
                size = 3, fontface = "bold", color = "black")
    
    print(p)
    
    ggsave(file.path(output_dir, paste0("feature_plot_", gene, "_with_labels.png")),
           p, width = 14, height = 8, dpi = 300)
  }
  
  dev.off()
  
  return(genes_present)
}

create_pathway_feature_plots <- function(seurat_obj, reduction = umap_reduction) {
  
  all_as_genes <- unlist(genes_of_interest)
  as_genes_present <- all_as_genes[all_as_genes %in% rownames(seurat_obj)]
  
  if(length(as_genes_present) > 0) {
    pdf(file.path(pathway_dir, "AS_genes_feature_plots_all.pdf"), width = 14, height = 8)
    
    for(i in seq_along(as_genes_present)) {
      gene <- as_genes_present[i]
      pathway <- names(which(sapply(genes_of_interest, function(x) gene %in% x)))[1]
      
      cat(paste("Plotting", gene, "from", pathway, "pathway (", i, "/", length(as_genes_present), ")...\n"))
      
      expr_data <- FetchData(seurat_obj, vars = c(gene, "predicted.celltype.l2", "condition", 
                                                  paste0(reduction, "_1"), paste0(reduction, "_2")))
      
      colnames(expr_data)[1] <- "gene_expression"
      
      expr_data <- expr_data %>%
        arrange(gene_expression)
      
      p <- ggplot(expr_data, aes_string(x = paste0(reduction, "_1"), 
                                        y = paste0(reduction, "_2"))) +
        geom_point(aes(color = gene_expression), size = 0.5) +
        scale_color_gradient(low = "lightgrey", high = "red", name = "Expression") +
        facet_wrap(~condition) +
        theme_minimal() +
        theme(
          strip.text = element_text(size = 18, face = "bold"),
          legend.position = "right",
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 16, hjust = 0.5)
        ) +
        ggtitle(gene, subtitle = paste("Pathway:", pathway))
      
      label_data <- expr_data %>%
        group_by(predicted.celltype.l2, condition) %>%
        summarise(
          x = median(get(paste0(reduction, "_1"))),
          y = median(get(paste0(reduction, "_2"))),
          .groups = "drop"
        )
      
      p <- p + 
        geom_text(data = label_data, 
                  aes(x = x, y = y, label = predicted.celltype.l2),
                  size = 3, fontface = "bold", color = "black")
      
      print(p)
    }
    
    dev.off()
    
    for(pathway_name in names(genes_of_interest)) {
      pathway_genes <- genes_of_interest[[pathway_name]]
      pathway_genes_present <- pathway_genes[pathway_genes %in% rownames(seurat_obj)]
      
      if(length(pathway_genes_present) > 0) {
        pdf(file.path(pathway_dir, paste0("feature_plots_", pathway_name, ".pdf")), 
            width = 14, height = 8)
        
        for(gene in pathway_genes_present) {
          expr_data <- FetchData(seurat_obj, vars = c(gene, "predicted.celltype.l2", "condition", 
                                                      paste0(reduction, "_1"), paste0(reduction, "_2")))
          
          colnames(expr_data)[1] <- "gene_expression"
          
          expr_data <- expr_data %>%
            arrange(gene_expression)
          
          p <- ggplot(expr_data, aes_string(x = paste0(reduction, "_1"), 
                                            y = paste0(reduction, "_2"))) +
            geom_point(aes(color = gene_expression), size = 0.5) +
            scale_color_gradient(low = "lightgrey", high = "red", name = "Expression") +
            facet_wrap(~condition) +
            theme_minimal() +
            theme(
              strip.text = element_text(size = 18, face = "bold"),
              legend.position = "right",
              axis.text = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              panel.grid = element_blank(),
              plot.title = element_text(size = 24, face = "bold", hjust = 0.5)
            ) +
            ggtitle(gene)
          
          label_data <- expr_data %>%
            group_by(predicted.celltype.l2, condition) %>%
            summarise(
              x = median(get(paste0(reduction, "_1"))),
              y = median(get(paste0(reduction, "_2"))),
              .groups = "drop"
            )
          
          p <- p + 
            geom_text(data = label_data, 
                      aes(x = x, y = y, label = predicted.celltype.l2),
                      size = 3, fontface = "bold", color = "black")
          
          print(p)
        }
        
        dev.off()
      }
    }
  }
}

# Violin plots
create_violin_plots <- function(seurat_obj, genes, group_by = "predicted.celltype.l2", 
                                split_by = "condition", test_method = "wilcox") {
  
  genes_present <- genes[genes %in% rownames(seurat_obj)]
  violin_plots <- list()
  stats_results <- data.frame()
  
  for(gene in genes_present) {
    p <- VlnPlot(seurat_obj, 
                 features = gene,
                 group.by = group_by,
                 split.by = split_by,
                 split.plot = TRUE,
                 cols = pub_colors$condition,
                 pt.size = 0,
                 adjust = 1) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        legend.position = "right",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16, face = "bold")
      ) +
      labs(title = gene, y = "Expression Level")
    
    violin_plots[[gene]] <- p
    
    celltypes <- unique(seurat_obj@meta.data[[group_by]])
    
    for(ct in celltypes) {
      cells_ct <- subset(seurat_obj, subset = predicted.celltype.l2 == ct)
      
      if(length(unique(cells_ct$condition)) == 2 && ncol(cells_ct) > 20) {
        expr_data <- LayerData(cells_ct, layer = "data", assay = "RNA")
        
        if(gene %in% rownames(expr_data)) {
          expr_as <- as.numeric(expr_data[gene, cells_ct$condition == "AS"])
          expr_hc <- as.numeric(expr_data[gene, cells_ct$condition == "Healthy"])
          
          if(length(expr_as) > 3 && length(expr_hc) > 3) {
            test_result <- wilcox.test(expr_as, expr_hc)
            
            stats_results <- rbind(stats_results, data.frame(
              gene = gene,
              celltype = ct,
              p_value = test_result$p.value,
              mean_AS = mean(expr_as),
              mean_HC = mean(expr_hc),
              log2FC = log2((mean(expr_as) + 0.1) / (mean(expr_hc) + 0.1)),
              n_AS = sum(cells_ct$condition == "AS"),
              n_HC = sum(cells_ct$condition == "Healthy")
            ))
          }
        }
      }
    }
  }
  
  if(nrow(stats_results) > 0) {
    stats_results$p_adj <- p.adjust(stats_results$p_value, method = "BH")
    stats_results$significance <- case_when(
      stats_results$p_adj < 0.001 ~ "***",
      stats_results$p_adj < 0.01 ~ "**",
      stats_results$p_adj < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  }
  
  saveRDS(list(plots = violin_plots, stats = stats_results), 
          file.path(plot_objects_dir, "violin_plots.rds"))
  
  return(list(plots = violin_plots, stats = stats_results))
}

create_gene_dotplot <- function(seurat_obj, genes, group_by = "condition") {
  genes_present <- genes[genes %in% rownames(seurat_obj)]
  
  if(length(genes_present) == 0) return(NULL)
  
  p_base <- DotPlot(seurat_obj, 
                    features = genes_present,
                    group.by = "predicted.celltype.l2",
                    split.by = group_by,
                    cols = c("lightblue", "darkred"),
                    dot.scale = 10,
                    scale = TRUE) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = "bold"),
      axis.text.y = element_text(size = 14),
      axis.title = element_blank(),
      legend.position = "right",
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.margin = margin(10, 10, 10, 50)
    ) +
    ggtitle("Gene Expression by Cell Type") +
    coord_flip()
  
  saveRDS(p_base, file.path(plot_objects_dir, "dotplot.rds"))
  
  return(p_base)
}

########################
######################## MAIN ANALYSIS

# DE analysis on all cells
de_all <- perform_de_analysis(integrated, output_prefix = "all_cells")
gsea_all <- perform_gsea(de_all, output_prefix = "all_cells")

# AS pathways
integrated <- analyze_as_pathways(integrated, de_all, output_prefix = "all_cells")

# cell type specific analyses
cd4_cells <- subset(integrated, subset = predicted.celltype.l2 %in% c("CD4 TCM", "CD4 TEM", "CD4 Naive"))
cat("\nCD4 T cells subset size:", ncol(cd4_cells), "\n")
if(ncol(cd4_cells) > 100) {
  de_cd4 <- perform_de_analysis(cd4_cells, output_prefix = "CD4_T_cells")
  gsea_cd4 <- perform_gsea(de_cd4, output_prefix = "CD4_T_cells")
  cd4_cells <- analyze_as_pathways(cd4_cells, de_cd4, output_prefix = "CD4_T_cells")
}

mono_cells <- subset(integrated, subset = predicted.celltype.l2 %in% c("CD14 Mono", "CD16 Mono"))
cat("\nMonocytes subset size:", ncol(mono_cells), "\n")
if(ncol(mono_cells) > 100) {
  de_mono <- perform_de_analysis(mono_cells, output_prefix = "Monocytes")
  gsea_mono <- perform_gsea(de_mono, output_prefix = "Monocytes")
  mono_cells <- analyze_as_pathways(mono_cells, de_mono, output_prefix = "Monocytes")
}

# feature plots
genes_to_plot <- key_genes[key_genes %in% rownames(integrated)]
cat("Genes found for plotting:", paste(genes_to_plot, collapse = ", "), "\n")
if(length(genes_to_plot) > 0) {
  create_individual_feature_plots(integrated, genes = genes_to_plot)
}
create_pathway_feature_plots(integrated)

# violin plots
violin_results <- create_violin_plots(integrated, genes = key_genes, test_method = "wilcox")

if(length(violin_results$plots) > 0) {
  pdf(file.path(output_dir, "violin_plots_key_genes.pdf"), width = 16, height = 10)
  for(gene in names(violin_results$plots)) {
    print(violin_results$plots[[gene]])
  }
  dev.off()
  
  for(gene in names(violin_results$plots)) {
    ggsave(file.path(output_dir, paste0("violin_plot_", gene, ".png")),
           violin_results$plots[[gene]], width = 16, height = 10, dpi = 300)
  }
}

if(nrow(violin_results$stats) > 0) {
  write.csv(violin_results$stats %>% 
              arrange(p_adj) %>%
              mutate(across(where(is.numeric), ~round(., 4))), 
            file.path(output_dir, "gene_expression_stats_by_celltype.csv"),
            row.names = FALSE)
  
  cat("\nTop differentially expressed genes by cell type:\n")
  top_findings <- violin_results$stats %>%
    filter(p_adj < 0.05) %>%
    arrange(p_adj) %>%
    head(10)
  print(top_findings)
}

# dotplot
dotplot_genes <- create_gene_dotplot(integrated, genes = key_genes)
if(!is.null(dotplot_genes)) {
  ggsave(file.path(output_dir, "gene_expression_dotplot.pdf"),
         dotplot_genes, width = 14, height = 16, dpi = 300)
  ggsave(file.path(output_dir, "gene_expression_dotplot.png"),
         dotplot_genes, width = 14, height = 16, dpi = 300)
}

# Mod score correlation
module_scores <- integrated@meta.data %>%
  dplyr::select(contains("_Score1")) %>%
  as.matrix()

if(ncol(module_scores) > 1) {
  cor_matrix <- cor(module_scores, use = "complete.obs")
  
  pdf(file.path(pathway_dir, "pathway_correlation_heatmap.pdf"), width = 10, height = 10)
  corrplot::corrplot(cor_matrix, 
                     method = "color",
                     type = "upper",
                     order = "hclust",
                     tl.col = "black",
                     tl.srt = 45,
                     addCoef.col = "black",
                     number.cex = 0.7,
                     col = pub_colors$heatmap)
  dev.off()
}

# summary figures
summary_stats <- integrated@meta.data %>%
  group_by(condition, predicted.celltype.l2) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(condition) %>%
  mutate(percent = n / sum(n) * 100)

p_summary <- ggplot(summary_stats, aes(x = predicted.celltype.l2, y = percent, fill = condition)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8) +
  scale_fill_manual(values = pub_colors$condition) +
  theme_cowplot() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    legend.position = "top",
    legend.text = element_text(size = 14),
    legend.title = element_blank(),
    plot.margin = margin(10, 10, 10, 50)
  ) +
  labs(x = "Cell Type", y = "Percentage of Cells", 
       title = "Cell Type Distribution by Condition")

ggsave(file.path(output_dir, "celltype_distribution.pdf"),
       p_summary, width = 14, height = 10, dpi = 300)
ggsave(file.path(output_dir, "celltype_distribution.png"),
       p_summary, width = 14, height = 10, dpi = 300)

if(exists("violin_results") && nrow(violin_results$stats) > 0) {
  sig_genes_summary <- violin_results$stats %>%
    filter(p_adj < 0.05) %>%
    group_by(gene) %>%
    summarise(
      n_sig_celltypes = n(),
      mean_log2FC = mean(log2FC),
      .groups = "drop"
    ) %>%
    arrange(desc(n_sig_celltypes), desc(abs(mean_log2FC)))
  
  if(nrow(sig_genes_summary) > 0) {
    p_sig_summary <- ggplot(sig_genes_summary %>% head(20), 
                            aes(x = reorder(gene, n_sig_celltypes), y = n_sig_celltypes)) +
      geom_bar(stat = "identity", fill = "#D62728") +
      coord_flip() +
      theme_cowplot() +
      theme(
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = "bold")
      ) +
      labs(x = "Gene", y = "Number of Cell Types with Significant Change",
           title = "Genes with Significant Changes Across Cell Types")
    
    ggsave(file.path(output_dir, "significant_genes_summary.pdf"),
           p_sig_summary, width = 12, height = 10, dpi = 300)
    ggsave(file.path(output_dir, "significant_genes_summary.png"),
           p_sig_summary, width = 12, height = 10, dpi = 300)
  }
}

# Combine pathway summaries across analyses
all_pathway_summaries <- data.frame()

for(analysis in c("all_cells", "CD4_T_cells", "Monocytes")) {
  summary_file <- file.path(pathway_dir, paste0("pathway_summary_comprehensive_", analysis, ".csv"))
  if(file.exists(summary_file)) {
    temp_summary <- read.csv(summary_file)
    temp_summary$analysis <- analysis
    all_pathway_summaries <- rbind(all_pathway_summaries, temp_summary)
  }
}

if(nrow(all_pathway_summaries) > 0) {
  p_compare <- all_pathway_summaries %>%
    ggplot(aes(x = pathway_display, y = percent_sig, fill = analysis)) +
    geom_bar(stat = "identity", position = "dodge") +
    coord_flip() +
    scale_fill_manual(values = c("all_cells" = "#1F77B4", 
                                 "CD4_T_cells" = "#FF7F0E", 
                                 "Monocytes" = "#2CA02C")) +
    theme_cowplot() +
    labs(x = "Pathway", y = "Percentage of Genes DE", 
         title = "AS Pathway Changes Across Cell Populations",
         fill = "Analysis") +
    theme(legend.position = "top")
  
  ggsave(file.path(pathway_dir, "pathway_comparison_across_analyses.pdf"),
         p_compare, width = 12, height = 8)
  
  write.csv(all_pathway_summaries,
            file.path(pathway_dir, "pathway_summary_all_analyses_combined.csv"),
            row.names = FALSE)
}

cat("Results saved to:", output_dir, "\n")
