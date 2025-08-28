#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(tidyverse)
    library(Seurat)
    library(SeuratData)
    library(Azimuth)
    library(SeuratDisk)
})

################################################################################
#   Utility Functions
################################################################################
# azimuth_ref_info <- function(refname) {
#     suppressPackageStartupMessages({
#         require(SeuratData)
#     })
#     # Load reference
#     SeuratData::InstallData(refname)
#     # reference <- SeuratData::LoadData(refname, type = "azimuth")$map
#     reference <- SeuratData::LoadData(refname)$map
#     annotation.levels <- names(slot(object = reference, name = "meta.data"))
#     annotation.levels <- annotation.levels[!grepl(pattern = "^nCount", x = annotation.levels)]
#     annotation.levels <- annotation.levels[!grepl(pattern = "^nFeature", x = annotation.levels)]
#     annotation.levels <- annotation.levels[!grepl(pattern = "^ori", x = annotation.levels)]
#     has.adt <- "ADT" %in% names(reference@assays)
#     return(list(annotation.levels, has.adt))
# }

my.AddAzimuthResults <- function(
        the.seurat, 
        the.seurat.mapped, 
        lvls, 
        tmpdir = tempdir()
) {
    # Adding Azimuth results using Seurat::AddAzimuthResults requires
    # an Azimuth mapping scores file, which is an RDS file containing a list
    # with:
    #   results$umap : Dimensionality reduction (UMAP) as a `DimReduc` object
    #   results$pred.df : `data.frame` with predicted celltypes and scores
    #   results$impADT :  `Assay` object with imputed values for antibody derived
    #                      tags (optional)
    results <- list()

    if('impADT' %in% Seurat::Assays(the.seurat.mapped)) {
        results$impADT <- the.seurat.mapped[['impADT']]
    }

    if('ref.umap' %in% Seurat::Reductions(the.seurat.mapped)) {
        results$umap <- the.seurat.mapped[['ref.umap']]
    }

    results$pred.df <- the.seurat.mapped@meta.data %>%
        tibble::rownames_to_column('cell') %>%
        dplyr::select(
            cell,
            dplyr::starts_with('predicted.'),
            mapping.score
        ) %>%
        as.data.frame

    az.outfile <- file.path(tmpdir, paste0('tmpAzimuthResults.rds'))
    saveRDS(results, file = az.outfile)

    # Add Azimuth results
    the.seurat <- Seurat::AddAzimuthResults(the.seurat, filename=az.outfile)
    file.remove(az.outfile)

    # Predicted celltype, score, and mapping.score columns are added to the
    # metadata but are all NA (bug?). Add it manually here
    stopifnot(all(rownames(the.seurat@meta.data) == results$pred.df$cell))
    for (n in lvls) {
        pcol <- paste0('predicted.', n)
        scol <- paste0('predicted.', n, '.score')
        the.seurat@meta.data[[pcol]] <- results$pred.df[[pcol]]
        the.seurat@meta.data[[scol]] <- results$pred.df[[scol]]
    }
    the.seurat@meta.data[['mapping.score']] <- results$pred.df[['mapping.score']]

    # Update assays
    for (n in lvls) {
        aname <- paste0('prediction.score.', n)
        the.seurat[[aname]] <- the.seurat.mapped@assays[[aname]]
    }

    # Return Seurat object with mapping information
    return(the.seurat)
}

################################################################################
#   Main Function
################################################################################

# annotate_azimuth.run <- function(
#         the.seurat, 
#         refname
# ) {
#     suppressPackageStartupMessages({
#         require(Azimuth)
#     })
#     # Get reference info
#     # refinfo <- azimuth_ref_info(refname)
#     # ann.levels <- refinfo[[1]]
#     # do.adt <- refinfo[[2]]
#     
#     # Run Azimuth
#     message("Running Azimuth...")
#     
#     ### This was needed because of expired https cert
#     # options(url.method='libcurl')
#     # httr::set_config(httr::config(ssl_verifypeer = 0L))
#     do.adt <- TRUE
#     print(mapped <- Azimuth::RunAzimuth(
#         the.seurat,
#         reference = refname,
#         do.adt = do.adt
#     ))
#     
#     # Add Azimuth results to Seurat object
#     my.AddAzimuthResults(the.seurat, mapped, ann.levels)
# }

annotate_unsupervised.run <- function(
        the.seurat, 
        ndim = 30,
        res = 0.8
) {
    # Run unsupervised workflow
    cat("Running unsupervised workflow...\n")
    
    the.seurat %>% 
        SCTransform() %>% 
        RunPCA(verbose=FALSE) %>% 
        RunUMAP(dims=1:ndim, verbose=FALSE) %>%
        FindNeighbors(dims=1:ndim, verbose=FALSE) %>%
        FindClusters(resolution = res, verbose=FALSE)
}



assign_idents.run <- function(
        counts.mtx,
        barcodes.tsv, 
        features.tsv,
        out.idents.tsv,
        out.inames.tsv,
        refname = NULL,
        local.ref = NULL,
        out.h5 = NULL
) {
    if(is.null(refname)) {
        message('Unsupervised clustering')
        q()
    } else {
        if(is.null(local.ref)) {
            message('Loading reference data for ', refname)
            SeuratData::InstallData(refname)
            refobj <- SeuratData::LoadData(refname)
            refarg <- refname
        } else {
            message('Loading local reference data for ', refname)
            refobj <- readRDS(file.path(local.ref, refname, 'ref.Rds'))
            refarg <- file.path(local.ref, refname)
        }
        do.adt <- "ADT" %in% names(refobj@assays)
        annotation.levels <- names(slot(object = refobj, name = "meta.data")) %>%
            stringr::str_subset("^nCount", negate=T) %>%
            stringr::str_subset("^nFeature", negate=T) %>%
            stringr::str_subset("^ori", negate=T)
        
        message("Loading count matrix...")
        print(sobj <- Seurat::CreateSeuratObject(
            Seurat::ReadMtx(
                mtx = counts.mtx,
                cells = barcodes.tsv,
                features = features.tsv
            )
        ))
        
        message("Running Azimuth...")
        print(mapped <- Azimuth::RunAzimuth(
            sobj,
            reference = refarg,
            do.adt = do.adt
        ))
        
        #--- Extract identities
        idents.df <- mapped@meta.data %>%
            select(starts_with('predicted') & !ends_with('score'))
        
        #--- Rename identities
        idents.names <- data.frame(
            row.names = names(idents.df),
            alt.name = str_c('l', 1:ncol(idents.df))
        )
        
        #--- Create output dir
        dir.create(dirname(out.idents.tsv), recursive=TRUE, showWarnings=FALSE)
        
        #--- Save idents as TSV
        write.table(
            idents.df,
            file=out.idents.tsv,
            sep = '\t', quote=F
        )
        
        #--- Save idents names
        write.table(
            idents.names, 
            file=out.inames.tsv,
            sep='\t', col.names=F, quote=F
        )
        
        #--- Create Seurat object with mapping data
        if(!is.null(out.h5)) {
            my.AddAzimuthResults(sobj, mapped, annotation.levels) %>%
                SeuratDisk::SaveH5Seurat(
                    filename=out.h5,
                    overwrite=T
                )
        }
    }
}

################################################################################
#   Parse arguments
################################################################################
if (sys.nframe() == 0L) {
    if(exists("snakemake")) {
        message('Seurat version ', packageVersion('Seurat'))
        message('Azimuth version ', packageVersion('Azimuth'))
        message('SeuratData version ', packageVersion('SeuratData'))
        assign_idents.run(
            counts.mtx = snakemake@input[['counts_mtx']],
            barcodes.tsv = snakemake@input[['barcodes_tsv']],
            features.tsv = snakemake@input[['features_tsv']],
            out.idents.tsv = snakemake@output[['out_idents_tsv']],
            out.inames.tsv = snakemake@output[['out_inames_tsv']],
            refname = snakemake@params[['refname']],
            local.ref = snakemake@params[['local_ref']],
            out.h5 = snakemake@output[['out_h5']]
        )
    } else {
        # args = commandArgs(trailingOnly=TRUE)
        # usage <- 'USAGE:\n  canonical_idents.R [-h] indir refname outdir'
        # if(args[1] == '-h' | args[1] == '--help') {
        #     cat(usage, '\n')
        #     q("no")
        # }
        # if(length(args) != 3)
        #     stop(usage, call.=FALSE)
        # 
        # outdir <- args[3]
        # if(!dir.exists(outdir)) 
        #     dir.create(outdir, recursive=TRUE, showWarnings=FALSE)        
        # 
        # assign_idents.run(
        #     counts.mtx = file.path(args[1], 'matrix.mtx'),
        #     barcodes.tsv = file.path(args[1], 'barcodes.tsv'),
        #     features.tsv = file.path(args[1], 'features.tsv'),
        #     refname = args[2],
        #     outdir = outdir
        # )
    }
}

cat("#--- assign_idents.R completed. ---#\n")


### testing
# suppressPackageStartupMessages({
#     require(Seurat)
#     require(dplyr)
#     require(magrittr)
# })
# 
# samp <- 'd4.U'
# counts.mtx <- sprintf('results/canonicalQC/%s/final/matrix.mtx', samp)
# barcodes.tsv <- sprintf('results/canonicalQC/%s/final/barcodes.tsv', samp)
# features.tsv <- sprintf('results/canonicalQC/%s/final/features.tsv', samp)
# refname <- 'pbmcref'
# # 
# print(sobj <- Seurat::CreateSeuratObject(
#     Seurat::ReadMtx(
#         mtx = counts.mtx,
#         cells = barcodes.tsv,
#         features = features.tsv
#     )
# ))
# 
# avail.datasets <- SeuratData::AvailableData()$Dataset
# refname %in% avail.datasets
# 
# sobj <- annotate_azimuth.run(sobj, refname)
# idents.df <- sobj@meta.data %>%
#     select(starts_with('predicted') & !ends_with('score'))
# 
# 
# 
# # barcodes.tsv <- sprintf('results/align_multi_starsolo/%s/Solo.out/Gene/raw/barcodes.tsv', samp)
# # features.tsv <- sprintf('results/align_multi_starsolo/%s/Solo.out/Gene/raw/features.tsv', samp)
