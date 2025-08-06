#!/usr/bin/env Rscript

# Load required libraries
library(Seurat)
library(Azimuth)
library(tidyverse)
library(patchwork)
library(future)
library(Matrix)
library(glmGamPoi)

options(timeout=400)

plan("multisession", workers = 4)

output_dir <- "results/b_integration_azimuth"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

all_samples <- readLines("config/all_samples_list")
hc_samples <- all_samples[grep("^HC", all_samples)]
as_samples <- all_samples[grep("^AS", all_samples)]

samples_to_analyze <- c(hc_samples, as_samples)

cat("Processing", length(hc_samples), "HC samples and", length(as_samples), "AS samples\n")
cat("Total samples to analyze:", length(samples_to_analyze), "\n\n")

read_sample_data <- function(sample_name, base_path = "results/b_demultiplexed_final") {
  sample_path <- file.path(base_path, sample_name)
  
  if (!dir.exists(sample_path)) {
    warning(paste("Directory not found:", sample_path))
    return(NULL)
  }
  
  mat <- ReadMtx(
    mtx = file.path(sample_path, "matrix.mtx"),
    features = file.path(sample_path, "features.tsv"),
    cells = file.path(sample_path, "barcodes.tsv"),
    feature.column = 1  
  )
  
  colnames(mat) <- paste0(sample_name, "_", colnames(mat))
  
  return(mat)
}

seurat_list <- list()

for (i in seq_along(samples_to_analyze)) {
  sample_name <- samples_to_analyze[i]
  cat("Processing", i, "/", length(samples_to_analyze), ":", sample_name, "\n")
  
  mat <- read_sample_data(sample_name)
  
  if (is.null(mat)) {
    cat("  Skipping due to read error\n")
    next
  }
  
  seurat_obj <- CreateSeuratObject(
    counts = mat,
    project = sample_name
  )
  
  # metadata
  seurat_obj$sample_id <- sample_name
  seurat_obj$condition <- ifelse(grepl("^HC", sample_name), "Healthy", "AS")
  seurat_obj$orig.barcode <- gsub(paste0("^", sample_name, "_"), "", colnames(seurat_obj))
  
  # Calculate mitochondrial percentage
  #seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, assay = "RNA", pattern = "^MT-")
  
  cat("  Cells:", ncol(seurat_obj), "| Features:", nrow(seurat_obj), "\n")
  cat("  Cell barcode example:", colnames(seurat_obj)[1], "\n")
  
  seurat_list[[sample_name]] <- seurat_obj
}

seurat_list <- seurat_list[!sapply(seurat_list, is.null)]

cat("\nSuccessfully loaded", length(seurat_list), "samples\n")

all_cell_names <- unlist(lapply(seurat_list, colnames))
if (length(all_cell_names) != length(unique(all_cell_names))) {
  stop("ERROR: Duplicate cell names found! Check barcode naming.")
} else {
  cat("✓ All cell names are unique\n")
}

################################ Annotate
# august 2025 This will only work if the Seurat version is 5.0.2, the joys of r dependencies
# I could cry right now
# test <- seurat_list$HC05
# test <- SCTransform(test)
# test <- Azimuth::RunAzimuth(test, reference="pbmcref")

process_seurat_list <- function(seurat_list, reference = "pbmcref", verbose = TRUE) {
  processed_list <- list()
  
  # Get names of objects
  object_names <- names(seurat_list)
  
  for (i in seq_along(seurat_list)) {
    obj_name <- object_names[i]
    
    if (verbose) {
      message(paste0("Processing ", obj_name, " (", i, "/", length(seurat_list), ")..."))
    }
    
    tryCatch({
      current_obj <- seurat_list[[i]]
      
      if (verbose) message("  Running SCTransform...")
      current_obj <- SCTransform(current_obj, verbose = FALSE)
      
      if (verbose) message("  Running Azimuth annotation...")
      current_obj <- Azimuth::RunAzimuth(current_obj, reference = reference)
      
      processed_list[[obj_name]] <- current_obj
      
      if (verbose) message(paste0("  ✓ ", obj_name, " completed successfully\n"))
      
    }, error = function(e) {
      warning(paste0("Error processing ", obj_name, ": ", e$message))
      processed_list[[obj_name]] <- NULL
    })
  }
  
  return(processed_list)
}

processed_seurat_list <- process_seurat_list(seurat_list, reference = "pbmcref")

dir.create("results/b_azimuth", showWarnings = FALSE, recursive = TRUE)

saveRDS(processed_seurat_list, file = "results/b_azimuth/seurat_list_filtered_azimuth.rds")
# note that three samples had less than 100 cells (cannot be mapped with azimuth)

source("scripts/plot_n_cells.R")





################################ Seurat v5 Integration Workflow
for(i in seq_along(processed_seurat_list)) {
  # Calculate from RNA assay
  processed_seurat_list[[i]][["percent.mt"]] <- PercentageFeatureSet(
    processed_seurat_list[[i]], 
    pattern = "^MT-",
    assay = "RNA"
  )
}

integrated <- merge(processed_seurat_list[[1]], 
                    y = processed_seurat_list[-1], 
                    add.cell.ids = names(processed_seurat_list),
                    project = "AS_HC_Integration")

# Verify metadata is preserved
print(paste("Number of metadata columns:", ncol(integrated@meta.data)))
print(colnames(integrated@meta.data))

# Since objects are already SCTransformed, we need to re-run SCTransform on the merged object
# First, set default assay to RNA
DefaultAssay(integrated) <- "RNA"

integrated <- SCTransform(integrated, 
                          vars.to.regress = "percent.mt", 
                          verbose = TRUE,
                          vst.flavor = "v2")

saveRDS(integrated, file.path(output_dir, "ckpt_integrated_sct_pre_azimuth.rds"))

# Run PCA
integrated <- RunPCA(integrated, npcs = 50, verbose = FALSE)

integrated <- IntegrateLayers(
  object = integrated,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony",
  verbose = TRUE
)

integrated[["RNA"]] <- JoinLayers(integrated[["RNA"]])

integrated <- FindNeighbors(integrated, reduction = "harmony", dims = 1:30)

plan("sequential")  
integrated <- FindClusters(integrated, resolution = 0.6, cluster.name = "harmony_clusters")

integrated <- RunUMAP(integrated, reduction = "harmony", dims = 1:30)

DefaultAssay(integrated) <- "RNA"
integrated <- NormalizeData(integrated, verbose = FALSE)

# metadata preservation
cat("\nFinal metadata check:\n")
print(paste("Total cells:", ncol(integrated)))
print(paste("Metadata columns:", ncol(integrated@meta.data)))
print(table(integrated$condition))
print(table(integrated$sample_id))

pred_cols <- grep("predicted", colnames(integrated@meta.data), value = TRUE)
cat("\nPrediction columns preserved:\n")
print(pred_cols)

saveRDS(integrated, file.path(output_dir, "integrated_seurat_harmony.rds"))
