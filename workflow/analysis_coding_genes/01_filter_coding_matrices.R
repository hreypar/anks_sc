# call libraries
library(Seurat)
#library(dplyr)
library(magrittr)
library(tidyverse)
require(DropletUtils)
library(scater)
library(Matrix)

# all Locally run

# functions
#
'%!in%' <- function(x,y)!('%in%'(x,y))

# set up the samples
c("PBMC-01-1-GEX-S1_SAMN25241888", "PBMC-01-2-GEX-S1_SAMN25241889", "PBMC-01-3-GEX-S1_SAMN25241890",
  "PBMC-01-4-GEX-S1_SAMN25241891", "PBMC-02-1-GEX-S2_SAMN25241930", "PBMC-02-2-GEX-S2_SAMN25241936",
  "PBMC-02-3-GEX-S2_SAMN25241934", "PBMC-02-4-GEX-S2_SAMN25241932", "PBMC-03-1-GEX-S3_SAMN25241924",
  "PBMC-03-2-GEX-S3_SAMN25241926", "PBMC-03-3-GEX-S3_SAMN25241928", "PBMC-03-4-GEX-S3_SAMN25241922",
  "PBMC-04-1-GEX-S4_SAMN25241905", "PBMC-04-2-GEX-S4_SAMN25241908", "PBMC-04-3-GEX-S4_SAMN25241912",
  "PBMC-04-4-GEX-S4_SAMN25241910", "PBMC-05-1-GEX-S5_SAMN25241918", "PBMC-05-2-GEX-S5_SAMN25241916",
  "PBMC-05-3-GEX-S5_SAMN25241920", "PBMC-05-4-GEX-S5_SAMN25241914", "PBMC-06-1-GEX-S7_SAMN25241898",
  "PBMC-06-2-GEX-S7_SAMN25241902", "PBMC-06-3-GEX-S7_SAMN25241900", "PBMC-06-4-GEX-S7_SAMN25241904",
  "PBMC-07-1-GEX-S7_SAMN25241892", "PBMC-07-2-GEX-S7_SAMN25241893", "PBMC-07-3-GEX-S7_SAMN25241897",
  "PBMC-07-4-GEX-S7_SAMN25241894",
  "AS1_SAMN20423166", "AS2_SAMN20423167", "AS3_SAMN20423168",
  "HC1_SAMN20423163", "HC2_SAMN20423164", "HC3_SAMN20423165") -> samples

# remove doublets from the liao samples (use the resources raw barcodes lists)
## read in liao exclude barcodes
## exclude from each sample
# empty drops
# is outlier
# write out seurat objects

# it has to be one by one bc I wont have enough memory to have them all at the same time

# Define directories
base_dir <- "results/align_multi_starsolo"
output_dir <- "results/b_canonicalQC"
barcodes_dir <- "config/barcodes_raw"


# Function to process each sample
process_sample <- function(sample, base_dir, output_dir, barcodes_dir,
#                           emptyDrops.niter=100000, emptyDrops.fdr.th = 1e-3,
                           min.upper.percent.mt = 5) {
 
   message(paste0("Processing sample: ", sample))
  
  # make output directory for the sample
  sample_output_dir <- file.path(output_dir, sample)
  dir.create(sample_output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Read input files  
#  sobj <- Seurat::ReadSTARsolo(data.dir = file.path(base_dir, sample, "Solo.out/Gene/raw"))
  sobj <- Seurat::ReadSTARsolo(data.dir = file.path(base_dir, sample, "Solo.out/Gene/filtered"))
  
  sobj <- Seurat::CreateSeuratObject(counts = sobj, project = sample)
  
  # genetic Doublet Removal (only for PBMC samples)
  if (grepl("^PBMC", sample)) {
    exclude_barcodes_file <- file.path(barcodes_dir, paste0("exclude_barcodes.", sample, ".tsv"))
    doublet_barcodes <- readLines(exclude_barcodes_file)
    keep_barcodes <- colnames(sobj) %!in% doublet_barcodes
    sobj <- sobj[, keep_barcodes]
  }
  
  # # empty Drops Removal using emptyDrops
  # ed.out <- DropletUtils::emptyDrops(
  #   sobj@assays$RNA@layers$counts,
  #   niters = emptyDrops.niter
  #   )
  # 
  # colnames(ed.out) <- str_c('edrop.', colnames(ed.out))
  # sobj@meta.data <- sobj@meta.data %>% cbind(ed.out)
  # 
  # sobj <- subset(sobj, subset = !(edrop.FDR < emptyDrops.fdr.th & edrop.Limited))
  
  # mitochondrial filtering and outlier detection using scater
  sobj[["percent.mt"]] <- Seurat::PercentageFeatureSet(sobj, pattern = "^MT-")
  
  sobj@meta.data %<>%
    dplyr::mutate(
      out.nCount_RNA = scater::isOutlier(nCount_RNA, log = TRUE, type = 'both'),
      out.nFeature_RNA = scater::isOutlier(nFeature_RNA, log = TRUE, type = 'both')
    )
  
  out.percent.mt <- scater::isOutlier(sobj$percent.mt, type = 'higher')
  upper.percent.mt <- max(min.upper.percent.mt, attr(out.percent.mt, 'thresholds')['higher'])
  
  # combine all outlier information 
  sobj@meta.data <- sobj@meta.data %>%
    dplyr::mutate(out.percent.mt = is.nan(percent.mt) | (percent.mt > upper.percent.mt)) %>%
    dplyr::mutate(discard = out.nCount_RNA | out.nFeature_RNA | out.percent.mt) # Combined outlier flags
  
  sobj <- subset(sobj, subset = !discard)
  
  # outputs
  # Save as rds
  rds_path <- file.path(sample_output_dir, paste0(sample, "_filtered.rds"))
  saveRDS(sobj, file = rds_path)
  
  # save as text-based sparse matrix for Scrublet
  mtx_filtered_path <- file.path(sample_output_dir, "matrix.mtx")
  barcodes_filtered_path <- file.path(sample_output_dir, "barcodes.tsv")
  features_filtered_path <- file.path(sample_output_dir, "features.tsv")
  
  Matrix::writeMM(sobj@assays$RNA$counts, file = mtx_filtered_path)
  write.table(Seurat::Cells(sobj), file = barcodes_filtered_path, quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(data.frame(row.names=rownames(sobj), gene_names=rownames(sobj)), file = features_filtered_path, quote = FALSE, row.names = FALSE, col.names=FALSE, sep="\t")
  
  
  # Remove objects to free memory
  rm(sobj)
  gc()
  
  message(paste0("Finished processing sample: ", sample))
}


# run them all 
for (sample in samples) {
  process_sample(sample, base_dir, output_dir, barcodes_dir,
#                 emptyDrops.niter = 10000,  
#                 emptyDrops.fdr.th = 0.01,
                 min.upper.percent.mt = 8)
} 
