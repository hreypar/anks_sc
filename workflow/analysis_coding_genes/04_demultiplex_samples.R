#!/usr/bin/env Rscript

# 04_demultiplex_samples.R
# Demultiplex pooled samples and organize all samples in a consistent structure

demultiplex_and_organize <- function(
    base_dir = "results/b_canonicalQC",
    metadata_file = "config/pbmc_anks_metadata.csv",
    barcodes_metadata_dir = "config/barcodes_metadata",
    output_dir = "results/b_demultiplexed_final",  # Changed from c_ to b_
    use_symlinks = TRUE
) {
  suppressPackageStartupMessages({
    require(Matrix)
    require(dplyr)
    require(tidyr)
    require(ggplot2)
  })
  
  # Read sample metadata
  sample_metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)
  cat(sprintf("Loaded metadata for %d samples\n", nrow(sample_metadata)))
  
  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Initialize summary
  summary_list <- list()
  
  # Process each sample
  for (i in 1:nrow(sample_metadata)) {
    sample_info <- sample_metadata[i, ]
    sample_name <- sample_info$sample
    sample_type <- sample_info$sample_type
    
    cat(sprintf("\n=== Processing %s (type: %s) ===\n", sample_name, sample_type))
    
    # Check if final filtered data exists
    input_dir <- file.path(base_dir, sample_name, "final")
    if (!dir.exists(input_dir)) {
      cat(sprintf("Warning: No final data found for %s, skipping...\n", sample_name))
      next
    }
    
    if (sample_type == "pool") {
      # Process pooled sample
      summary <- process_pooled_sample(
        sample_name = sample_name,
        input_dir = input_dir,
        barcodes_metadata_dir = barcodes_metadata_dir,
        output_dir = output_dir,
        base_dir = base_dir
      )
      summary_list <- c(summary_list, summary)
      
    } else if (sample_type == "sample") {
      # Process single sample (create symlink or copy)
      summary <- process_single_sample(
        sample_name = sample_name,
        input_dir = input_dir,
        output_dir = output_dir,
        use_symlinks = use_symlinks
      )
      summary_list[[length(summary_list) + 1]] <- summary
      
    } else {
      cat(sprintf("Unknown sample type: %s\n", sample_type))
    }
  }
  
  # Create summary report and plots
  create_summary_report(summary_list, output_dir)
  
  cat("\n=== Demultiplexing complete! ===\n")
  cat(sprintf("Results saved to: %s\n", output_dir))
}

process_pooled_sample <- function(
    sample_name,
    input_dir,
    barcodes_metadata_dir,
    output_dir,
    base_dir
) {
  cat("Processing pooled sample...\n")
  
  # Read barcode metadata
  metadata_file <- file.path(barcodes_metadata_dir, 
                             sprintf("barcodes_metadata.%s.csv", sample_name))
  
  if (!file.exists(metadata_file)) {
    stop(sprintf("Barcode metadata not found: %s", metadata_file))
  }
  
  barcode_metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)
  cat(sprintf("  Loaded metadata for %d barcodes\n", nrow(barcode_metadata)))
  
  # Load matrix data
  cat("  Loading matrix data...\n")
  counts_mtx <- file.path(input_dir, "matrix.mtx")
  barcodes_tsv <- file.path(input_dir, "barcodes.tsv")
  features_tsv <- file.path(input_dir, "features.tsv")
  
  spmat <- Matrix::readMM(counts_mtx)
  barcodes <- readLines(barcodes_tsv)
  features_df <- read.table(features_tsv, header = FALSE, sep = '\t', 
                            stringsAsFactors = FALSE, fill = TRUE)
  
  colnames(spmat) <- barcodes
  rownames(spmat) <- make.unique(features_df$V1)
  
  cat(sprintf("  Loaded matrix: %d features x %d cells\n", nrow(spmat), ncol(spmat)))
  
  # Match barcodes
  matched_metadata <- barcode_metadata %>%
    filter(barcode %in% barcodes)
  
  cat(sprintf("  Matched %d barcodes between matrix and metadata\n", nrow(matched_metadata)))
  
  # Check for unmatched barcodes in matrix
  unmatched_barcodes <- barcodes[!barcodes %in% barcode_metadata$barcode]
  if (length(unmatched_barcodes) > 0) {
    cat(sprintf("  WARNING: %d barcodes in matrix not found in metadata!\n", 
                length(unmatched_barcodes)))
    
    # Save unmatched barcodes to file
    unmatched_file <- file.path(output_dir, 
                                sprintf("unmatched_barcodes_%s.txt", sample_name))
    writeLines(unmatched_barcodes, unmatched_file)
    cat(sprintf("  Unmatched barcodes saved to: %s\n", unmatched_file))
  }
  
  # Check for barcodes in metadata but not in matrix (informational)
  metadata_only_barcodes <- barcode_metadata$barcode[!barcode_metadata$barcode %in% barcodes]
  if (length(metadata_only_barcodes) > 0) {
    cat(sprintf("  INFO: %d barcodes in metadata not found in matrix (likely filtered earlier)\n", 
                length(metadata_only_barcodes)))
  }
  
  # Get unique subjects
  subjects <- unique(matched_metadata$Subject)
  cat(sprintf("  Found %d unique subjects\n", length(subjects)))
  
  # Initialize summary
  summary_list <- list()
  
  # Add summary for unmatched barcodes if any
  if (length(unmatched_barcodes) > 0) {
    summary_list[[length(summary_list) + 1]] <- list(
      sample = paste0(sample_name, "_UNMATCHED"),
      original_pool = sample_name,
      n_cells = length(unmatched_barcodes),
      n_features = nrow(spmat),
      type = "unmatched",
      status = "no_metadata"
    )
  }
  
  # Split by subject
  for (subject in subjects) {
    cat(sprintf("\n  Processing subject: %s\n", subject))
    
    # Get barcodes for this subject
    subject_barcodes <- matched_metadata %>%
      filter(Subject == subject) %>%
      pull(barcode)
    
    cat(sprintf("    %d cells for this subject\n", length(subject_barcodes)))
    
    # Subset matrix
    subject_mat <- spmat[, subject_barcodes, drop = FALSE]
    
    # Create output directory for subject
    subject_dir <- file.path(output_dir, subject)
    dir.create(subject_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Write output files
    cat("    Writing output files...\n")
    
    # Matrix
    Matrix::writeMM(subject_mat, file = file.path(subject_dir, "matrix.mtx"))
    
    # Barcodes
    writeLines(colnames(subject_mat), file.path(subject_dir, "barcodes.tsv"))
    
    # Features (same for all subjects)
    write.table(features_df, file = file.path(subject_dir, "features.tsv"),
                sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    # Add metadata file with original pool information
    subject_metadata <- matched_metadata %>%
      filter(Subject == subject) %>%
      mutate(original_pool = sample_name)
    
    write.csv(subject_metadata, 
              file = file.path(subject_dir, "cell_metadata.csv"),
              row.names = FALSE)
    
    # Summary
    summary_list[[length(summary_list) + 1]] <- list(
      sample = subject,
      original_pool = sample_name,
      n_cells = length(subject_barcodes),
      n_features = nrow(subject_mat),
      type = "demultiplexed",
      status = matched_metadata %>% 
        filter(Subject == subject) %>% 
        pull(Status) %>% 
        unique() %>% 
        paste(collapse = ",")
    )
    
    cat(sprintf("    ✓ Saved to: %s\n", subject_dir))
  }
  
  # Create reconciliation report
  recon_file <- file.path(output_dir, sprintf("reconciliation_%s.txt", sample_name))
  separator <- paste(rep("=", nchar(sample_name) + 25), collapse = "")
  recon_lines <- c(
    sprintf("Reconciliation Report for %s", sample_name),
    separator,
    sprintf(""),
    sprintf("Matrix barcodes: %d", length(barcodes)),
    sprintf("Metadata barcodes: %d", nrow(barcode_metadata)),
    sprintf("Matched barcodes: %d", nrow(matched_metadata)),
    sprintf("Unmatched in matrix: %d", length(unmatched_barcodes)),
    sprintf("Metadata-only barcodes: %d", length(metadata_only_barcodes)),
    sprintf(""),
    sprintf("Coverage: %.1f%% of matrix barcodes found in metadata", 
            100 * nrow(matched_metadata) / length(barcodes))
  )
  writeLines(recon_lines, recon_file)
  
  # Create barplot for this pool - SAVE TO FINAL DIRECTORY
  create_pool_barplot(matched_metadata, sample_name, input_dir)  # Changed from base_dir to input_dir
  
  return(summary_list)
}

create_pool_barplot <- function(matched_metadata, sample_name, final_dir) {
  cat("  Creating cell distribution barplot...\n")
  
  # Summarize cells per subject
  subject_summary <- matched_metadata %>%
    group_by(Subject, Status) %>%
    summarise(n_cells = n(), .groups = "drop") %>%
    arrange(desc(n_cells))
  
  # Create barplot
  p <- ggplot(subject_summary, aes(x = reorder(Subject, n_cells), y = n_cells, fill = Status)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = n_cells), hjust = -0.2) +
    coord_flip() +
    labs(title = paste("Cell Distribution -", sample_name),
         subtitle = paste("Total cells:", sum(subject_summary$n_cells)),
         x = "Subject",
         y = "Number of Cells") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.text = element_text(size = 10)
    ) +
    scale_fill_brewer(palette = "Set2")
  
  # Save to final directory (which is passed as final_dir parameter)
  plot_file <- file.path(final_dir, "demultiplexing_distribution.pdf")
  ggsave(plot_file, p, width = 8, height = 6)
  cat(sprintf("  Barplot saved to: %s\n", plot_file))
}

process_single_sample <- function(
    sample_name,
    input_dir,
    output_dir,
    use_symlinks = TRUE
) {
  cat("Processing single sample...\n")
  
  # For single samples, use the sample name as is
  sample_dir <- file.path(output_dir, sample_name)
  dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Get absolute paths for symlinks
  input_dir_abs <- normalizePath(input_dir)
  
  # Files to link/copy
  files <- c("matrix.mtx", "barcodes.tsv", "features.tsv")
  
  for (file in files) {
    input_file <- file.path(input_dir_abs, file)
    output_file <- file.path(sample_dir, file)
    
    if (use_symlinks) {
      # Create symlink
      if (!file.exists(output_file)) {
        file.symlink(input_file, output_file)
        cat(sprintf("  Created symlink: %s\n", file))
      }
    } else {
      # Copy file
      file.copy(input_file, output_file, overwrite = TRUE)
      cat(sprintf("  Copied: %s\n", file))
    }
  }
  
  # Count cells
  n_cells <- length(readLines(file.path(sample_dir, "barcodes.tsv")))
  
  # Create metadata file for consistency
  metadata_df <- data.frame(
    sample = sample_name,
    type = "single_sample",
    n_cells = n_cells
  )
  write.csv(metadata_df, 
            file = file.path(sample_dir, "sample_info.csv"),
            row.names = FALSE)
  
  cat(sprintf("  ✓ Saved to: %s\n", sample_dir))
  
  return(list(
    sample = sample_name,
    original_pool = NA,
    n_cells = n_cells,
    n_features = NA,  # Could count if needed
    type = "single_sample",
    status = "original"
  ))
}

create_summary_report <- function(summary_list, output_dir) {
  cat("\nCreating summary report...\n")
  
  # Convert to data frame
  summary_df <- do.call(rbind, lapply(summary_list, as.data.frame))
  
  # Save detailed summary
  summary_file <- file.path(output_dir, "demultiplexing_summary.tsv")
  write.table(summary_df, file = summary_file, sep = "\t", 
              quote = FALSE, row.names = FALSE)
  
  # Print summary statistics
  cat("\n=== Summary Statistics ===\n")
  cat(sprintf("Total entries: %d\n", nrow(summary_df)))
  cat(sprintf("  - Demultiplexed subjects: %d\n", 
              sum(summary_df$type == "demultiplexed")))
  cat(sprintf("  - Single samples: %d\n", 
              sum(summary_df$type == "single_sample")))
  cat(sprintf("  - Unmatched barcode groups: %d\n", 
              sum(summary_df$type == "unmatched", na.rm = TRUE)))
  cat(sprintf("Total cells: %d\n", sum(summary_df$n_cells)))
  
  # Check for unmatched barcodes
  if (any(summary_df$type == "unmatched", na.rm = TRUE)) {
    cat("\n⚠️  WARNING: Some barcodes could not be matched to metadata!\n")
    unmatched_summary <- summary_df %>%
      filter(type == "unmatched") %>%
      select(original_pool, n_cells)
    print(unmatched_summary)
    
    cat("\nCheck the following files for details:\n")
    unmatched_files <- list.files(output_dir, pattern = "unmatched_barcodes.*\\.txt", 
                                  full.names = TRUE)
    cat(paste("  -", unmatched_files, collapse = "\n"), "\n")
  }
  
  # Summary by pool
  if (any(summary_df$type == "demultiplexed")) {
    cat("\nDemultiplexing summary by pool:\n")
    pool_summary <- summary_df %>%
      filter(type == "demultiplexed") %>%
      group_by(original_pool) %>%
      summarise(
        n_subjects = n(),
        total_cells = sum(n_cells),
        .groups = "drop"
      )
    print(pool_summary)
    
    # Create faceted plot for all pools
    create_all_pools_plot(summary_df, output_dir)
  }
  
  cat(sprintf("\nDetailed summary saved to: %s\n", summary_file))
}

create_all_pools_plot <- function(summary_df, output_dir) {
  cat("\nCreating combined pools visualization...\n")
  
  # Filter for demultiplexed samples only
  demux_data <- summary_df %>%
    filter(type == "demultiplexed") %>%
    mutate(
      pool_short = gsub("_SAMN.*", "", original_pool),
      Subject = sample
    )
  
  # Create faceted barplot
  p_facet <- ggplot(demux_data, aes(x = reorder(Subject, n_cells), y = n_cells, fill = status)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = n_cells), hjust = -0.2, size = 3) +
    coord_flip() +
    facet_wrap(~ pool_short, scales = "free", ncol = 2) +
    labs(title = "Cell Distribution Across All Pooled Samples",
         subtitle = paste("Total demultiplexed cells:", sum(demux_data$n_cells)),
         x = "Subject",
         y = "Number of Cells",
         fill = "Status") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12),
      strip.text = element_text(size = 10, face = "bold"),
      axis.text = element_text(size = 8),
      legend.position = "bottom"
    ) +
    scale_fill_brewer(palette = "Set2")
  
  # Save faceted plot
  facet_file <- file.path(output_dir, "all_pools_distribution.pdf")
  ggsave(facet_file, p_facet, width = 14, height = 10)
  cat(sprintf("Faceted plot saved to: %s\n", facet_file))
  
  # Create summary barplot by pool
  pool_summary <- demux_data %>%
    group_by(pool_short) %>%
    summarise(
      total_cells = sum(n_cells),
      n_subjects = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(total_cells))
  
  p_summary <- ggplot(pool_summary, aes(x = reorder(pool_short, total_cells), y = total_cells)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_text(aes(label = paste0(total_cells, "\n(", n_subjects, " subjects)")), 
              hjust = -0.1, size = 3.5) +
    coord_flip() +
    labs(title = "Total Cells per Pool",
         subtitle = "After QC and demultiplexing",
         x = "Pool",
         y = "Total Cells") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 10)
    )
  
  summary_file <- file.path(output_dir, "pools_summary.pdf")
  ggsave(summary_file, p_summary, width = 8, height = 6)
  cat(sprintf("Summary plot saved to: %s\n", summary_file))
}

# Command line interface
if (sys.nframe() == 0L) {
  suppressPackageStartupMessages({
    require(optparse)
  })
  
  option_list <- list(
    make_option(c("-b", "--base-dir"), type = "character", 
                default = "results/b_canonicalQC",
                help = "Base directory with filtered samples [default: %default]"),
    make_option(c("-m", "--metadata"), type = "character", 
                default = "config/pbmc_anks_metadata.csv",
                help = "Sample metadata file [default: %default]"),
    make_option(c("-d", "--barcode-dir"), type = "character", 
                default = "config/barcodes_metadata",
                help = "Directory with barcode metadata files [default: %default]"),
    make_option(c("-o", "--output-dir"), type = "character", 
                default = "results/b_demultiplexed_final",  # Changed from c_ to b_
                help = "Output directory [default: %default]"),
    make_option(c("--copy"), action = "store_true", default = FALSE,
                help = "Copy files instead of creating symlinks for single samples")
  )
  
  parser <- OptionParser(option_list = option_list,
                         description = "Demultiplex pooled samples and organize all samples")
  args <- parse_args(parser)
  
  # Run demultiplexing
  demultiplex_and_organize(
    base_dir = args$`base-dir`,
    metadata_file = args$metadata,
    barcodes_metadata_dir = args$`barcode-dir`,
    output_dir = args$`output-dir`,
    use_symlinks = !args$copy
  )
}
