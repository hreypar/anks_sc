#!/usr/bin/env Rscript

filter.run <- function(
        counts.mtx,
        barcodes.tsv,
        features.tsv,
        include.barcodes.tsv,
        out.counts.mtx,
        out.barcodes.tsv,
        out.features.tsv,
        out.excluded.tsv = NULL,
        out.raw.h5 = NULL,
        out.filtered.h5 = NULL,
        barcodes.column = 1,
        emptyDrops.niter = 100000,
        emptyDrops.fdr.th = 1e-3,
        min.upper.percent.mt = 5,
        BPPARAM = BiocParallel::SerialParam()
) {
    suppressPackageStartupMessages({
        require(magrittr)
        require(stringr)
        require(dplyr)
        require(Seurat)
        require(DropletUtils)
    })
    
    #--- Create raw Seurat object
    print(sobj <- Seurat::CreateSeuratObject(
        Seurat::ReadMtx(
            mtx = counts.mtx,
            cells = barcodes.tsv,
            features = features.tsv
        )
    ))
    
    #--- Add barcodes to include from external file
    if(file.exists(include.barcodes.tsv)) {
        include.barcodes <- read.table(
            include.barcodes.tsv, header=FALSE, sep='\t'
        )
        sobj@meta.data$ext.include <- 
            Seurat::Cells(sobj) %in% include.barcodes[, barcodes.column]
        message(sprintf('ext.include = %d cells', sum(sobj$ext.include)))
    } else {
        sobj@meta.data$ext.include <-FALSE
    }
    
    #--- EmptyDrops
    ed.out <- DropletUtils::emptyDrops(
        sobj[['RNA']]@layers$counts,
        niters=emptyDrops.niter,
        BPPARAM = BPPARAM
    )
    colnames(ed.out) <- str_c('edrop.', colnames(ed.out))
    sobj@meta.data %<>% 
        cbind(ed.out) %>%
        dplyr::mutate(
            edrop.is_cell = !is.na(edrop.FDR) & edrop.FDR < emptyDrops.fdr.th
        )
    message(sprintf('EmptyDrops identified %d cells', sum(sobj$edrop.is_cell)))
    
    #--- Add percent mitochondrial
    sobj@meta.data %<>%
        dplyr::mutate(
            percent.mt = Seurat::PercentageFeatureSet(sobj, pattern = '^MT-')
        )
    
    #--- Detect outliers for nCount and nFeature
    sobj@meta.data %<>%
        dplyr::mutate(
            out.nCount_RNA = scater::isOutlier(
                nCount_RNA, 
                subset=edrop.is_cell,
                log=TRUE,
                type='both'
            ),
            out.nFeature_RNA = scater::isOutlier(
                nFeature_RNA,
                subset=edrop.is_cell,
                log=TRUE, 
                type='both'
            )
        )
    
    #--- Detect outliers for percent.mt
    out.percent.mt <- scater::isOutlier(
        sobj$percent.mt,
        subset = sobj$edrop.is_cell,
        type='higher'
    )
    
    upper.percent.mt <- max(min.upper.percent.mt, attr(out.percent.mt, 'thresholds')['higher'])
    message(sprintf('upper bound for percent.mt = %.2f', upper.percent.mt))
    sobj@meta.data %<>%
        dplyr::mutate(
            out.percent.mt = is.nan(percent.mt) | (percent.mt > upper.percent.mt)
        )
    
    #--- Final pass: "included" or an EmptyDrops cell that is not an outlier
    #    for nCount, nFeature, or percent.mt
    sobj@meta.data %<>%
        dplyr::mutate(
            pass.qc = !(out.nCount_RNA | out.nFeature_RNA | out.percent.mt)
        ) %>%
        dplyr::mutate(
            pass.final = ext.include | (edrop.is_cell & pass.qc)
        )
    
    message(sprintf('raw cells: %s', dim(sobj)[2]))
    message(sprintf('non-empty cells: %s', sum(sobj$edrop.is_cell)))
    message(sprintf('passing QC: %s cells', sum(sobj$pass.qc)))
    message(sprintf('non-empty passing QC: %s cells', sum(sobj$edrop.is_cell & sobj$pass.qc)))
    message(sprintf('pass.final: %d cells', sum(sobj$pass.final)))
    message(sprintf('included but not cell: %d cells', sum(sobj$ext.include & (!sobj$edrop.is_cell))))
    message(sprintf('included but not passing QC: %d cells', sum(sobj$ext.include & (!sobj$pass.qc))))  
    
    #--- Subset Seurat object
    print(sobj.final <- sobj[,sobj$pass.final])

    #--- Create output dir
    dir.create(dirname(out.counts.mtx), recursive=TRUE, showWarnings=FALSE)
    
    #--- Output MTX
    Matrix::writeMM(sobj.final[['RNA']]@layers$counts, file=out.counts.mtx)
    
    #--- Output barcodes
    Seurat::Cells(sobj.final) %>%
        write.table(
            file = out.barcodes.tsv,
            sep='\t', quote=F, row.names=F, col.names=F
        )
    
    #--- Output excluded barcodes
    if(!is.null(out.excluded.tsv)) {
        setdiff(Seurat::Cells(sobj), Seurat::Cells(sobj.final)) %>%
            write.table(
                file = out.excluded.tsv,
                sep='\t', quote=F, row.names=F, col.names=F
            )
    }
    
    #--- Output features
    features.df <- read.table(features.tsv, header=F, sep='\t')
    stopifnot(all(Features(sobj) == make.unique(features.df$V2)))
    features.df %>%
        write.table(
            file = out.features.tsv,
            sep='\t', quote=F, row.names=F, col.names=F
        )
    
    if(!is.null(out.raw.h5)) {
        message('writing raw h5 file to ', out.raw.h5)
        SeuratDisk::SaveH5Seurat(
            sobj,
            filename = out.raw.h5, 
            overwrite = T
        )
    }
    
    if(!is.null(out.filtered.h5)) {
        message('writing filtered h5 file to ', out.filtered.h5)
        SeuratDisk::SaveH5Seurat(
            sobj.final, 
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
        if(snakemake@threads > 1) {
            bpparam <- BiocParallel::MulticoreParam(workers=snakemake@threads)
        } else {
            bpparam <- BiocParallel::SerialParam()
        }
        filter.run(
            snakemake@input[["counts_mtx"]],
            snakemake@input[["barcodes_tsv"]],
            snakemake@input[["features_tsv"]],
            snakemake@input[["include_barcodes_tsv"]],
            snakemake@output[["counts_mtx"]],
            snakemake@output[["barcodes_tsv"]],
            snakemake@output[["features_tsv"]],
            out.excluded.tsv = snakemake@output[['excluded_tsv']],
            out.raw.h5 = snakemake@output[['raw_h5']],
            out.filtered.h5 = snakemake@output[['filtered_h5']],
            barcodes.column = snakemake@params[['barcodes_column']],
            emptyDrops.niter = snakemake@params[['emptyDrops_niter']],
            emptyDrops.fdr.th = snakemake@params[['emptyDrops_fdr_th']],
            min.upper.percent.mt = snakemake@params[['min_upper_percent_mt']],
            BPPARAM = bpparam
        )
    } else {
        samp <- 'd4.U'
        filter.run(
            sprintf('results/align_multi_starsolo/%s/Solo.out/Gene/raw/matrix.mtx', samp),
            sprintf('results/align_multi_starsolo/%s/Solo.out/Gene/raw/barcodes.tsv', samp),
            sprintf('results/align_multi_starsolo/%s/Solo.out/Gene/raw/features.tsv', samp),
            sprintf('resources/barcodes_raw/include_barcodes.%s.tsv', samp),
            sprintf('results/canonicalQC/%s/filtered/matrix.mtx'),
            sprintf('results/canonicalQC/%s/filtered/barcodes.tsv'),
            sprintf('results/canonicalQC/%s/filtered/features.tsv'),            
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
    message("#--- filter_raw.R completed. ---#\n")
}

################################################################################
#   testing
################################################################################
# samp <- 'd4.U'
# counts.mtx <- sprintf('results/align_multi_starsolo/%s/Solo.out/Gene/raw/matrix.mtx', samp)
# barcodes.tsv <- sprintf('results/align_multi_starsolo/%s/Solo.out/Gene/raw/barcodes.tsv', samp)
# features.tsv <- sprintf('results/align_multi_starsolo/%s/Solo.out/Gene/raw/features.tsv', samp)
# include.barcodes.tsv <- sprintf('resources/barcodes_raw/include_barcodes.%s.tsv', samp)
# out.raw.h5 <- sprintf('results/canonicalQC/%s/filtered/raw.h5seurat', samp)
# out.filtered.h5 <- sprintf('results/canonicalQC/%s/filtered/filtered.h5seurat', samp)
# out.counts.mtx <- sprintf('results/canonicalQC/%s/filtered/matrix.mtx', samp)
# out.barcodes.tsv <- sprintf('results/canonicalQC/%s/filtered/barcodes.tsv', samp)
# out.features.tsv <- sprintf('results/canonicalQC/%s/filtered/features.tsv', samp)

