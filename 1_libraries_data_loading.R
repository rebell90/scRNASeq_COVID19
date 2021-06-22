

####### Load Packages -- includes automatic installation from CRAN/Bioconductor if package is not yet installed
#install.packages("pacman")
pacman::p_load(scRNAseq, Glimma, DropletUtils,iSEE,BiocParallel,mclust,dplyr,biomaRt, scater,robustbase, scran, BiocNeighbors,
               Seurat, patchwork,htmltools, raster, readr, data.table, stringr, ggplot2, SingleCellExperiment, sctransform, 
               AnnotationHub, limma, edgeR, tidyverse, cowplot, Matrix.utils, magrittr, Matrix, purrr, reshape2, S4Vectors,
               tibble, pheatmap, apeglm, png, DESeq2, ashr, sva, PCAtools, GEOquery, RColorBrewer) 

gc() ## garbage cleanup
setwd('/Users/beccabell/RNASeq/https:/github.com/rebell90/covid_rna_seq_project.git') ## set wd
mem.maxVSize(vsize = 100000000000) ##vector memory size

## download file from GEO GSE161918(listed as first supp-file)
covid_seurat <- readRDS( file = "GSE161918_AllBatches_SeuratObj.rds")

####### Single Cell Experiment ######

# Extract raw counts and metadata to create SingleCellExperiment object.
# data will be restricted to only the first time-point in the first batch
counts <- covid_seurat@assays$RNA@counts ## RNA-seq count matrix 
cluster_id <- DataFrame(covid_seurat@active.ident) ## cell-type 
meta_data <- covid_seurat@meta.data[which(covid_seurat@meta.data$Timepoint == "T0" &  
                                            str_trim(covid_seurat@meta.data$Batch) == 'B1'), ]


## subset count data so only cells pulled from T0/B1are in matrix
counts <- counts[, which(colnames(counts) %in% rownames(meta_data))]

## extract row names of metadata to accurately subset cluster_id values
cell_labels <- rownames(meta_data)
cluster_id <- cluster_id[cell_labels,]
meta_data <- cbind(meta_data, cluster_id) ## join cluster-id to metadata object

## Delete original seurat-object to open memory/storage
# covid_seurat <- NULL


# Create single cell experiment object
sce.covid <- SingleCellExperiment(assays = list(counts = counts), 
                                  colData = meta_data)

## view data object
sce.covid
