#!/usr/bin/env Rscript

# 03_final_filtering.R
# Remove doublets identified by Scrublet and save filtered 

scrublet_filter <- function(
    input_dir,
    output_dir = NULL,
    feature_col = 1,  # column contains feature names
    generate_qc = TRUE
) {
  suppressPackageStartupMessages({
    require(Matrix)
    require(dplyr)
    require(ggplot2)
  })
  
  if (is.null(output_dir)) {
    output_dir <- file.path(input_dir, "final")
  }
  
  counts.mtx <- file.path(input_dir, "matrix.mtx")
  barcodes.tsv <- file.path(input_dir, "barcodes.tsv")
  features.tsv <- file.path(input_dir, "features.tsv")
  pass.barcodes.tsv <- file.path(input_dir, "scrublet", "SCR.barcodes.pass.tsv")
  predictions.all.tsv <- file.path(input_dir, "scrublet", "SCR.predictions_all.tsv")
  
  required_files <- c(counts.mtx, barcodes.tsv, features.tsv, pass.barcodes.tsv)
  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files) > 0) {
    stop("Missing required files:\n", paste("  -", missing_files, collapse = "\n"))
  }
  
  sample_name <- basename(input_dir)
  cat(sprintf("\n=== Processing %s ===\n", sample_name))
  
  cat("Loading data...\n")
  spmat <- Matrix::readMM(counts.mtx)
  barcodes <- readLines(barcodes.tsv)
  features.df <- read.table(features.tsv, header = FALSE, sep = '\t', 
                            stringsAsFactors = FALSE, fill = TRUE)
  
  colnames(spmat) <- barcodes
  
  if (feature_col > ncol(features.df)) {
    stop(sprintf("Feature column %d does not exist (file has %d columns)", 
                 feature_col, ncol(features.df)))
  }
  rownames(spmat) <- make.unique(features.df[[feature_col]])
  
  cat(sprintf("Original matrix: %d features x %d cells\n", nrow(spmat), ncol(spmat)))
  
  # passing barcodes
  pass.barcodes <- readLines(pass.barcodes.tsv)
  cat(sprintf("Passing barcodes: %d\n", length(pass.barcodes)))
  
  if (!all(pass.barcodes %in% colnames(spmat))) {
    missing <- sum(!pass.barcodes %in% colnames(spmat))
    stop(sprintf("Error: %d passing barcodes not found in matrix", missing))
  }
  
  n_original <- ncol(spmat)
  n_passing <- length(pass.barcodes)
  n_removed <- n_original - n_passing
  pct_removed <- round(100 * n_removed / n_original, 2)
  
  cat(sprintf("Removing %d cells (%.2f%%)\n", n_removed, pct_removed))
  
  if (generate_qc && file.exists(predictions.all.tsv)) {
    cat("Generating QC report...\n")
    
    predictions <- read.table(predictions.all.tsv, header = TRUE, 
                              sep = '\t', stringsAsFactors = FALSE)
    
    pdf_file <- file.path(input_dir, "scrublet", "filtering_qc.pdf")
    pdf(pdf_file, width = 10, height = 6)
    
    p1 <- ggplot(predictions, aes(x = scrublet_score, fill = scrublet_prediction)) +
      geom_histogram(alpha = 0.7, bins = 50, position = "identity") +
      scale_fill_manual(values = c("Singlet" = "#4CAF50", "Doublet" = "#F44336")) +
      labs(title = paste(sample_name, "- Scrublet Score Distribution"),
           subtitle = sprintf("Removed %d doublets (%.1f%%)", n_removed, pct_removed),
           x = "Doublet Score", y = "Count") +
      theme_minimal() +
      theme(legend.position = "top")
    print(p1)
    
    dev.off()
    cat(sprintf("QC report saved to: %s\n", pdf_file))
  }
  
  spmat <- spmat[, pass.barcodes]
  
  cat(sprintf("Final matrix: %d features x %d cells\n", nrow(spmat), ncol(spmat)))
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  cat("Writing filtered data...\n")
  
  out.counts.mtx <- file.path(output_dir, "matrix.mtx")
  out.barcodes.tsv <- file.path(output_dir, "barcodes.tsv")
  out.features.tsv <- file.path(output_dir, "features.tsv")
  
  Matrix::writeMM(spmat, file = out.counts.mtx)
  
  writeLines(colnames(spmat), out.barcodes.tsv)
  
  # preserve original structure (!)
  write.table(features.df, file = out.features.tsv,
              sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  summary_file <- file.path(output_dir, "filtering_summary.txt")
  # Fixed the separator line
  separator <- paste(rep("=", nchar(sample_name) + 30), collapse = "")
  summary_lines <- c(
    sprintf("Scrublet Filtering Summary for %s", sample_name),
    separator,
    sprintf("Date: %s", Sys.Date()),
    sprintf(""),
    sprintf("Original cells: %d", n_original),
    sprintf("Doublets removed: %d (%.2f%%)", n_removed, pct_removed),
    sprintf("Cells retained: %d", n_passing),
    sprintf(""),
    sprintf("Feature column used: %d", feature_col),
    sprintf("Input directory: %s", input_dir),
    sprintf("Output directory: %s", output_dir)
  )
  writeLines(summary_lines, summary_file)
  
  cat(sprintf("\nFiltering completed! Results saved to: %s\n", output_dir))
  
  invisible(list(
    sample = sample_name,
    n_original = n_original,
    n_removed = n_removed,
    n_retained = n_passing,
    pct_removed = pct_removed
  ))
}

if (sys.nframe() == 0L) {
  suppressPackageStartupMessages({
    require(optparse)
  })
  
  option_list <- list(
    make_option(c("-i", "--input-dir"), type = "character", 
                help = "Input directory containing matrix files and scrublet results"),
    make_option(c("-o", "--output-dir"), type = "character", default = NULL,
                help = "Output directory (default: input_dir/final)"),
    make_option(c("-f", "--feature-col"), type = "integer", default = 1,
                help = "Column number containing feature names (default: 1)"),
    make_option(c("--no-qc"), action = "store_true", default = FALSE,
                help = "Skip QC plot generation"),
    make_option(c("-a", "--all"), action = "store_true", default = FALSE,
                help = "Process all samples in results/b_canonicalQC/")
  )
  
  parser <- OptionParser(
    option_list = option_list,
    description = "Filter cells based on Scrublet doublet predictions",
    usage = "%prog [options]"
  )
  args <- parse_args(parser)
  
  if (args$all) {
    base_dir <- "results/b_canonicalQC"
    sample_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
    
    all_results <- list()
    
    for (sample_dir in sample_dirs) {
      # Check if scrublet results exist
      if (file.exists(file.path(sample_dir, "scrublet", "SCR.barcodes.pass.tsv"))) {
        tryCatch({
          result <- scrublet_filter(
            input_dir = sample_dir,
            feature_col = args$`feature-col`,
            generate_qc = !args$`no-qc`
          )
          all_results[[length(all_results) + 1]] <- result
        }, error = function(e) {
          cat(sprintf("Error processing %s: %s\n", basename(sample_dir), e$message))
        })
      } else {
        cat(sprintf("Skipping %s - no scrublet results found\n", basename(sample_dir)))
      }
    }
    
    if (length(all_results) > 0) {
      summary_df <- do.call(rbind, lapply(all_results, as.data.frame))
      summary_file <- file.path(base_dir, "scrublet_filtering_summary.tsv")
      write.table(summary_df, file = summary_file, sep = "\t", 
                  quote = FALSE, row.names = FALSE)
      cat(sprintf("\nSummary saved to: %s\n", summary_file))
      print(summary_df)
    }
    
  } else {
    # Process single sample
    if (is.null(args$`input-dir`)) {
      print_help(parser)
      stop("Input directory is required")
    }
    
    scrublet_filter(
      input_dir = args$`input-dir`,
      output_dir = args$`output-dir`,
      feature_col = args$`feature-col`,
      generate_qc = !args$`no-qc`
    )
  }
}
