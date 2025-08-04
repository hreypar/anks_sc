#!/bin/bash

# this script assumes
# the python script is at scripts/run_scrublet.py
# input files are named matrix.mtx, features.tsv, and barcodes.tsv
# we have the scrublet conda environment set up

# Run Scrublet on samples after genetic demultiplexing and qc filters

# Exit on error
set -e

# Script configuration
PYTHON_SCRIPT="scripts/run_scrublet.py"
BASE_DIR="results/b_canonicalQC"
CONDA_ENV="scrublet"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")

# Main log directory for summary only
SUMMARY_LOG_DIR="logs/scrublet"
mkdir -p "$SUMMARY_LOG_DIR"
SUMMARY_LOG="$SUMMARY_LOG_DIR/summary_${TIMESTAMP}.log"

# Doublet rates after genetic demultiplexing
PBMC_DOUBLET_RATE=0.03  # Conservative estimate for remaining doublets
OTHER_DOUBLET_RATE=0.08

# Function to log messages to summary log
log_summary() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] $1" | tee -a "$SUMMARY_LOG"
}

# Activate conda environment
#log_summary "Activating conda environment: $CONDA_ENV"
#conda activate $CONDA_ENV

# Check if python script exists
if [ ! -f "$PYTHON_SCRIPT" ]; then
    log_summary "Error: Python script not found at $PYTHON_SCRIPT"
    exit 1
fi

log_summary "Starting Scrublet analysis on genetically demultiplexed samples"
log_summary "=================================================="
log_summary "PBMC samples (post-genetic demux): doublet rate = $PBMC_DOUBLET_RATE"
log_summary "Other samples (single-sample): doublet rate = $OTHER_DOUBLET_RATE"
log_summary ""

# Track statistics
total_samples=0
processed_samples=0
skipped_samples=0
failed_samples=0

# Process each sample directory
for sample_dir in "$BASE_DIR"/*; do
    # Check if it's a directory
    if [ ! -d "$sample_dir" ]; then
        continue
    fi
    
    total_samples=$((total_samples + 1))
    
    # Get sample name
    sample_name=$(basename "$sample_dir")
    
    # Check if required input files exist
    if [ ! -f "$sample_dir/matrix.mtx" ] || \
       [ ! -f "$sample_dir/features.tsv" ] || \
       [ ! -f "$sample_dir/barcodes.tsv" ]; then
        log_summary "Warning: Missing input files for $sample_name, skipping..."
        skipped_samples=$((skipped_samples + 1))
        continue
    fi
    
    # Skip if already processed (remove this block to reprocess)
    if [ -f "$sample_dir/scrublet/SCR.summary.txt" ]; then
        log_summary "Sample $sample_name already processed, skipping..."
        skipped_samples=$((skipped_samples + 1))
        continue
    fi
    
    # Count current cells
    current_cells=$(wc -l < "$sample_dir/barcodes.tsv")
    
    # Determine doublet rate based on sample type
    if [[ $sample_name == PBMC* ]]; then
        doublet_rate=$PBMC_DOUBLET_RATE
        sample_type="PBMC (demultiplexed)"
    else
        doublet_rate=$OTHER_DOUBLET_RATE
        sample_type="AS/HC (single-sample)"
    fi
    
    log_summary "Processing $sample_name"
    log_summary "  Sample type: $sample_type"
    log_summary "  Current cells: $current_cells"
    log_summary "  Expected doublet rate: $doublet_rate"
    
    # Create output directory
    output_dir="$sample_dir/scrublet"
    mkdir -p "$output_dir"
    
    # Sample-specific log file in the scrublet directory
    sample_log="$output_dir/scrublet_run_${TIMESTAMP}.log"
    
    # Log sample info to the sample-specific log
    {
        echo "Scrublet Analysis Log"
        echo "===================="
        echo "Date: $(date)"
        echo "Sample: $sample_name"
        echo "Sample type: $sample_type"
        echo "Input cells: $current_cells"
        echo "Expected doublet rate: $doublet_rate"
        echo ""
        echo "Processing output:"
        echo "-------------------"
    } > "$sample_log"
    
    # Run scrublet and capture output to both sample log and console
    python "$PYTHON_SCRIPT" \
        --counts-matrix "$sample_dir/matrix.mtx" \
        --features "$sample_dir/features.tsv" \
        --barcodes "$sample_dir/barcodes.tsv" \
        --output-dir "$output_dir" \
        --doublet-rate "$doublet_rate" \
        2>&1 | tee -a "$sample_log"
    
    # Check if outputs were created and add summary to log
    if [ -f "$output_dir/SCR.barcodes.pass.tsv" ]; then
        log_summary "✓ Successfully processed $sample_name"
        processed_samples=$((processed_samples + 1))
        
        # Add summary to sample log
        {
            echo ""
            echo "-------------------"
            echo "Processing completed successfully"
            echo "Output files created:"
            ls -la "$output_dir"/SCR.* | awk '{print "  - " $9}'
            echo ""
            echo "Results summary:"
            cat "$output_dir/SCR.summary.txt"
        } >> "$sample_log"
        
        # Quick stats for summary
        detected_doublets=$(grep "Doublet" "$output_dir/SCR.summary.txt" | awk '{print $2}' || echo "0")
        log_summary "  Detected $detected_doublets potential doublets"
    else
        log_summary "✗ Error processing $sample_name"
        failed_samples=$((failed_samples + 1))
        echo "ERROR: Processing failed!" >> "$sample_log"
    fi
    
    echo "----------------------------------------"
done

# Generate summary report
log_summary "Generating summary report..."
summary_file="$SUMMARY_LOG_DIR/summary_results_${TIMESTAMP}.tsv"
{
    echo -e "Sample\tExperiment_Type\tCells_Input\tExpected_Rate\tDetected_Doublets\tDetected_Singlets\tDetected_Rate\tRemoval_Rate\tLog_File"
    for sample_dir in "$BASE_DIR"/*; do
        if [ -d "$sample_dir" ] && [ -f "$sample_dir/scrublet/SCR.summary.txt" ]; then
            sample_name=$(basename "$sample_dir")
            summary_file_sample="$sample_dir/scrublet/SCR.summary.txt"
            log_file="$sample_dir/scrublet/scrublet_run_${TIMESTAMP}.log"
            
            # Determine sample type
            if [[ $sample_name == PBMC* ]]; then
                sample_type="Multiplexed"
                expected_rate=$PBMC_DOUBLET_RATE
            else
                sample_type="Single-sample"
                expected_rate=$OTHER_DOUBLET_RATE
            fi
            
            # Get cell count and results
            if [ -f "$sample_dir/barcodes.tsv" ]; then
                current_cells=$(wc -l < "$sample_dir/barcodes.tsv")
                
                doublets=$(grep "Doublet" "$summary_file_sample" | awk '{print $2}' || echo "0")
                singlets=$(grep "Singlet" "$summary_file_sample" | awk '{print $2}' || echo "0")
                total=$((doublets + singlets))
                detected_rate=$(awk "BEGIN {printf \"%.4f\", $doublets/$total}" || echo "0")
                removal_rate=$(awk "BEGIN {printf \"%.1f%%\", ($doublets/$total)*100}" || echo "0%")
                
                # Relative path to log file
                rel_log_path="${sample_name}/scrublet/scrublet_run_${TIMESTAMP}.log"
                
                echo -e "$sample_name\t$sample_type\t$current_cells\t$expected_rate\t$doublets\t$singlets\t$detected_rate\t$removal_rate\t$rel_log_path"
            fi
        fi
    done
} > "$summary_file"

log_summary "Summary report saved to: $summary_file"

# Final summary
log_summary ""
log_summary "Processing Complete!"
log_summary "==================="
log_summary "Total samples found: $total_samples"
log_summary "Successfully processed: $processed_samples"
log_summary "Skipped: $skipped_samples"
log_summary "Failed: $failed_samples"

# Display results
echo ""
echo "Summary of results:"
column -t -s $'\t' "$summary_file" | head -20

# Create an index file listing all sample logs
index_file="$SUMMARY_LOG_DIR/sample_logs_index_${TIMESTAMP}.txt"
{
    echo "Sample Log Files Index"
    echo "====================="
    echo ""
    for sample_dir in "$BASE_DIR"/*; do
        if [ -f "$sample_dir/scrublet/scrublet_run_${TIMESTAMP}.log" ]; then
            sample_name=$(basename "$sample_dir")
            echo "$sample_name: $sample_dir/scrublet/scrublet_run_${TIMESTAMP}.log"
        fi
    done
} > "$index_file"

log_summary ""
log_summary "All logs saved:"
log_summary "  - Summary log: $SUMMARY_LOG"
log_summary "  - Results table: $summary_file"
log_summary "  - Sample logs index: $index_file"
log_summary "  - Individual sample logs: in each sample's scrublet/ directory"

