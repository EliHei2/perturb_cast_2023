# 2020-04-15 08:58
# elihei [<eheidari@student.ethz.ch>]
# /Volumes/Projects/scGCN/code/R/dep.R 

import_cran <- function(x) {
    if ( !require(x, character.only=TRUE) ) {
        install.packages(x, dep=TRUE, repos='http://cran.us.r-project.org')
        if ( !require(x, character.only=TRUE) ) 
            stop( "Package not found" )
    }
}

import_bioc <- function(x) {
    if ( !require(x, character.only=TRUE) ) {
        BiocManager::install(x)
        if ( !require(x, character.only=TRUE) ) 
            stop("Package not found")
    }
}

import_git  <- function(x) {
    if ( !require(x, character.only=TRUE) ) {
        devtools::install_git(x)
        if ( !require(x, character.only=TRUE) ) 
            stop( "Package not found" )
    }
}

import_int  <- function(x) {
    if ( !require(x, character.only=TRUE) ) {
        devtools::install('code/pika')
        if ( !require(x, character.only=TRUE) ) 
            stop( "Package not found" )
    }
}

library('magrittr')
import_cran('tidyverse')
import_cran('BiocManager')
import_cran('devtools')

cran_packages = c('gRim', 'glasso', 'methods', 'Matrix',
  'igraph', 'graph', 'data.table', 'visNetwork',
  'Rtsne', 'RColorBrewer', 'doParallel', 'reshape2', 'qgraph',
  'magrittr', 'uwot', 'scales', 'patchwork')

bioc_packages = c('genefilter', 'limma', 'BiocParallel',
  'gRbase', 'RBGL', 'Rgraphviz', 'gRain',
  'scran', 'Seurat', 'SpatialExperiment')

git_packages  = c('ropenscilabs/umapr')

int_packages  = c('pika')

cran_packages %>% map(import_cran)
bioc_packages %>% map(import_bioc)
# int_packages  %>% map(import_int)
# git_packages  %>% map(import_git)


# 09:48




