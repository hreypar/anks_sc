#! /usr/bin/env Rscript
print(.libPaths())

## Set default repo
local({
    r <- getOption("repos")
    r["CRAN"] <- "https://cloud.r-project.org"
    options(repos=r)
})

install.cran <- function(pkg, version=NULL) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        if(is.null(version)) {
            cat(sprintf('install.packages("%s", quiet=TRUE)\n', pkg))
            install.packages(pkg, quiet=TRUE)
        } else {
            cat(sprintf('remotes::install_version("%s", version="%s", upgrade=FALSE, quiet=TRUE)\n', pkg, version))
            remotes::install_version(pkg, version=version, upgrade=FALSE, quiet=TRUE)
        }
    } else {
        cat(sprintf('package "%s" already installed\n', pkg))
    }
}

install.bioc <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        cat(sprintf('BiocManager::install("%s", update=FALSE, quiet=TRUE)\n', pkg))
        BiocManager::install(pkg, update=FALSE, quiet=TRUE)
    } else {
        cat(sprintf('package "%s" already installed\n', pkg))
    }
}

install.github <- function(pkg, repo) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        cat(sprintf('remotes::install_github("%s", upgrade=FALSE, quiet=TRUE)\n', repo))
        remotes::install_github(repo, upgrade=FALSE, quiet=TRUE)
    } else {
        cat(sprintf('package "%s" already installed\n', pkg))
    }
}

install.cran('tidyverse')
install.cran('remotes')
install.cran('BiocManager')
install.cran('hdf5r')
install.cran('igraph')
install.bioc('SingleCellExperiment')
install.bioc("TFBSTools")
install.cran('SeuratObject')
install.cran('Seurat')
install.github('SeuratDisk', 'mojaveazure/seurat-disk')
install.github('SeuratData', 'satijalab/seurat-data')
install.github('Azimuth', 'satijalab/azimuth@master')
install.bioc('scater')
install.bioc('DropletUtils')

if(exists("snakemake")) {
    write.table(
        installed.packages()[,c(1,3:4)],
        file=snakemake@output[[1]],
        row.names=F, na = '', quote=F, sep='\t'
    )
}
