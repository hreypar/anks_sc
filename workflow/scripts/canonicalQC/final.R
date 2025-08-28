#!/usr/bin/env Rscript

final.run <- function(
        counts.mtx,
        barcodes.tsv,
        features.tsv,
        pass.barcodes.tsv,
        out.counts.mtx,
        out.barcodes.tsv,
        out.features.tsv,
        out.filtered.h5 = NULL
) {
    suppressPackageStartupMessages({
        require(magrittr)
        require(dplyr)
        require(Seurat)
    })
    
    #--- Create sparse matrix
    spmat <- Seurat::ReadMtx(
        mtx = counts.mtx,
        cells = barcodes.tsv,
        features = features.tsv
    )

    message('original cells: ', dim(spmat)[2])
    
    #--- Get passing barcodes
    barcodes.column <- 1
    pass.barcodes <- read.table(pass.barcodes.tsv,
                                header = FALSE, sep = '\t') %>%
        dplyr::pull(barcodes.column)
    message('passing barcodes: ', length(pass.barcodes))
    stopifnot(all(pass.barcodes %in% colnames(spmat)))
    
    #--- subset
    spmat <- spmat[ ,colnames(spmat) %in% pass.barcodes]
    
    message('final matrix: ', dim(spmat)[1], ' feats ', dim(spmat)[2], ' cells')
    
    #--- Create output dir
    dir.create(dirname(out.counts.mtx), recursive=TRUE, showWarnings=FALSE)
    
    #--- Output MTX
    Matrix::writeMM(spmat, file=out.counts.mtx)
    
    #--- Output barcodes
    colnames(spmat) %>%
        write.table(
            file = out.barcodes.tsv,
            sep='\t', quote=F, row.names=F, col.names=F
        )
    
    #--- Output features
    features.df <- read.table(features.tsv, header=F, sep='\t')
    stopifnot(all(rownames(spmat) == make.unique(features.df$V2)))
    features.df %>%
        write.table(
            file = out.features.tsv,
            sep='\t', quote=F, row.names=F, col.names=F
        )

    if(!is.null(out.filtered.h5)) {
        message('writing filtered h5 file to ', out.filtered.h5)
        Seurat::CreateSeuratObject(spmat) %>%
            SeuratDisk::SaveH5Seurat(
                filename = out.filtered.h5,
                overwrite = T
            )
    }
}


################################################################################
#   Parse arguments
################################################################################
if (sys.nframe() == 0L) {
    if(exists("snakemake")) {
        # if(snakemake@threads > 1) {
        #     bpparam <- BiocParallel::MulticoreParam(workers=snakemake@threads)
        # } else {
        #     bpparam <- BiocParallel::SerialParam()
        # }
        final.run(
            snakemake@input[["counts_mtx"]],
            snakemake@input[["barcodes_tsv"]],
            snakemake@input[["features_tsv"]],
            snakemake@input[["passBC_tsv"]],
            snakemake@output[["counts_mtx"]],
            snakemake@output[["barcodes_tsv"]],
            snakemake@output[["features_tsv"]],
            out.filtered.h5 = snakemake@output[['filtered_h5']]
        )
    } else {
        samp <- 'd4.U'
        final.run(
            sprintf('results/canonicalQC/%s/filtered/matrix.mtx', samp),
            sprintf('results/canonicalQC/%s/filtered/barcodes.tsv', samp),
            sprintf('results/canonicalQC/%s/filtered/features.tsv', samp),
            sprintf('results/canonicalQC/%s/scrublet/SCR.barcodes.pass.tsv', samp),
            sprintf('results/canonicalQC/%s/final/matrix.mtx', samp),
            sprintf('results/canonicalQC/%s/final/barcodes.tsv', samp),
            sprintf('results/canonicalQC/%s/final/features.tsv', samp)
        )
        
        # args = commandArgs(trailingOnly=TRUE)
        # usage <- 'USAGE:\n  cell_qc.R [-h] counts.dir outdir'
        # if(args[1] == '-h' | args[1] == '--help') {
        #     cat(usage, '\n')
        #     q("no")
        # }
        # if(length(args) != 2)
        #     stop(usage, call.=FALSE)
        # 
        # counts.mtx <- file.path(args[1], 'matrix.mtx')
        # barcodes.tsv <- file.path(args[1], 'barcodes.tsv')
        # features.tsv <- file.path(args[1], 'features.tsv')
        # stopifnot(file.exists(counts.mtx, barcodes.tsv, features.tsv))
        # outdir <- args[2]
    }
    message("#--- final.R completed. ---#\n")
}
