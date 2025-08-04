#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Scrublet wrapper for doublet detection in single-cell RNA-seq data.

This script runs Scrublet on filtered count matrices and outputs predictions
for doublet detection along with visualization plots.
"""

import argparse
import sys
import os
from pathlib import Path
from typing import Optional
from numbers import Number

import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def scrublet_wrapper(
        counts_fn: str,
        features_fn: str,
        barcodes_fn: str,
        outdir: str,
        fcol: int = 0,
        bcol: int = 0,
        expected_doublet_rate: Optional[float] = None, 
        sim_doublet_ratio: Number = 2.0,
        pass_barcodes_file: str = 'SCR.barcodes.pass.tsv',
        all_predictions_file: str = 'SCR.predictions_all.tsv',
        single_predictions_file: str = 'SCR.predictions_single.tsv',
        summary_file: str = 'SCR.summary.txt'
    ):
    """
    Run Scrublet doublet detection on single-cell RNA-seq data.
    
    Parameters
    ----------
    counts_fn : str
        Path to the counts matrix file (MTX format)
    features_fn : str
        Path to the features/genes file (TSV format)
    barcodes_fn : str
        Path to the barcodes file (TSV format)
    outdir : str
        Output directory for results
    fcol : int
        Column index for feature names (default: 0)
    bcol : int
        Column index for barcode names (default: 0)
    expected_doublet_rate : float, optional
        Expected doublet rate. If None, will be calculated based on cell count
    sim_doublet_ratio : float
        Ratio of simulated doublets to observed cells (default: 2.0)
    pass_barcodes_file : str
        Filename for passing barcodes output
    all_predictions_file : str
        Filename for all predictions output
    single_predictions_file : str
        Filename for singlet predictions output
    summary_file : str
        Filename for summary statistics output
    """
    
    # Load data
    print(f"Loading count matrix from {counts_fn}...", file=sys.stderr)
    counts_csc = scipy.io.mmread(counts_fn).T.tocsc()
    
    print(f"Loading features from {features_fn}...", file=sys.stderr)
    with open(features_fn, 'r') as fh:
        features = pd.Series(l.strip('\n').split('\t')[fcol] for l in fh)

    print(f"Loading barcodes from {barcodes_fn}...", file=sys.stderr)
    with open(barcodes_fn, 'r') as fh:
        barcodes = pd.Series(l.strip('\n').split('\t')[bcol] for l in fh)

    assert counts_csc.shape == (len(barcodes), len(features)), \
        f'Incongruent input: matrix shape {counts_csc.shape} != ({len(barcodes)}, {len(features)})'
    
    # Calculate expected doublet rate if not provided
    if expected_doublet_rate is None:
        # This is the calculation for 10x doublet rate but will be different for other platforms
        expected_doublet_rate = counts_csc.shape[0]/1000 * 0.008 

    print(f'Using expected doublet rate = {expected_doublet_rate}', file=sys.stderr)
    
    # Initialize Scrublet object
    scrub = scr.Scrublet(
        counts_csc,
        expected_doublet_rate=expected_doublet_rate,
        sim_doublet_ratio=sim_doublet_ratio,
    )

    # Run doublet detection
    print("Running doublet detection...", file=sys.stderr)
    _scores, _preds = scrub.scrub_doublets(
        min_counts=3,
        min_cells=3, 
        min_gene_variability_pctl=85, 
        n_prin_comps=30
    )

    # Create output directory
    os.makedirs(outdir, exist_ok=True)
    
    # Generate and save plots
    print("Generating plots...", file=sys.stderr)
    scrub.plot_histogram()
    plt.savefig(os.path.join(outdir, 'SCR.doublet_score_histogram.png'))
    plt.close()
    
    print('Running UMAP...', file=sys.stderr)
    scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
    print('Done.', file=sys.stderr)
    
    scrub.plot_embedding('UMAP', order_points=True)
    plt.savefig(os.path.join(outdir, 'SCR.UMAP.png'))
    plt.close()

    # Prepare results dataframe
    df = pd.DataFrame({
        'barcode': barcodes,
        'scrublet_score': pd.Series(scrub.doublet_scores_obs_),
        'pbool': pd.Series(scrub.predicted_doublets_)
    }).set_index('barcode')
    
    df['scrublet_prediction'] = 'Singlet'
    df.loc[df['pbool'], 'scrublet_prediction'] = 'Doublet'
    
    # Save outputs
    print("Saving results...", file=sys.stderr)
    
    # All predictions
    with open(os.path.join(outdir, all_predictions_file), 'w') as outh:
        df[['scrublet_prediction', 'scrublet_score']].to_csv(outh, sep='\t')
    
    # Singlet predictions only
    with open(os.path.join(outdir, single_predictions_file), 'w') as outh:
        df.loc[~df['pbool']][['scrublet_prediction', 'scrublet_score']].to_csv(outh, sep='\t')
    
    # Passing barcodes
    with open(os.path.join(outdir, pass_barcodes_file), 'w') as outh:
        pd.Series(df.loc[~df['pbool']].index).to_csv(outh, index=False, header=False)
    
    # Summary statistics
    summary = pd.DataFrame(df.scrublet_prediction.value_counts())
    with open(os.path.join(outdir, summary_file), 'w') as outh:
        summary.to_csv(outh, sep="\t")
    
    print(f"Results saved to {outdir}", file=sys.stderr)
    print(f"Total cells: {len(df)}", file=sys.stderr)
    print(f"Predicted doublets: {df['pbool'].sum()} ({df['pbool'].sum()/len(df)*100:.2f}%)", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description='Run Scrublet doublet detection on single-cell RNA-seq data',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Required arguments
    parser.add_argument('--counts-matrix', '-c', required=True,
                        help='Path to counts matrix file (MTX format)')
    parser.add_argument('--features', '-f', required=True,
                        help='Path to features/genes file (TSV format)')
    parser.add_argument('--barcodes', '-b', required=True,
                        help='Path to barcodes file (TSV format)')
    parser.add_argument('--output-dir', '-o', required=True,
                        help='Output directory for results')
    
    # Optional arguments
    parser.add_argument('--doublet-rate', '-r', type=float, default=None,
                        help='Expected doublet rate. If not provided, will be calculated automatically')
    parser.add_argument('--sim-doublet-ratio', type=float, default=2.0,
                        help='Ratio of simulated doublets to observed cells')
    parser.add_argument('--feature-column', type=int, default=0,
                        help='Column index for feature names in features file')
    parser.add_argument('--barcode-column', type=int, default=0,
                        help='Column index for barcode names in barcodes file')
    
    # Output file names
    parser.add_argument('--pass-barcodes-file', default='SCR.barcodes.pass.tsv',
                        help='Filename for passing barcodes output')
    parser.add_argument('--all-predictions-file', default='SCR.predictions_all.tsv',
                        help='Filename for all predictions output')
    parser.add_argument('--single-predictions-file', default='SCR.predictions_single.tsv',
                        help='Filename for singlet predictions output')
    parser.add_argument('--summary-file', default='SCR.summary.txt',
                        help='Filename for summary statistics output')
    
    args = parser.parse_args()
    
    # Run scrublet wrapper
    scrublet_wrapper(
        counts_fn=args.counts_matrix,
        features_fn=args.features,
        barcodes_fn=args.barcodes,
        outdir=args.output_dir,
        fcol=args.feature_column,
        bcol=args.barcode_column,
        expected_doublet_rate=args.doublet_rate,
        sim_doublet_ratio=args.sim_doublet_ratio,
        pass_barcodes_file=args.pass_barcodes_file,
        all_predictions_file=args.all_predictions_file,
        single_predictions_file=args.single_predictions_file,
        summary_file=args.summary_file
    )


if __name__ == '__main__':
    main()


