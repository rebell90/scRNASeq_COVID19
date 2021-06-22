

###### Quality Control #######
## With Single Cell data, it is very common to have sparse matrices, where there are many 0's.
## Most of these lowly expressed genes don't reveal much information , and create extra noise and take up space.
## There are an abundance of methods to filter lowly-expressed genes, and their selection is determined by analysis-goals.
## remove any genes that have a count of '0' in all cells

## First Step will be to filter out cells with spike-ins and a high-mitochondria percentage, as these indicate low-quality cells/reads
## Create Boolean Var to identify 
is.spike <- grepl("^ERCC", rownames(sce.covid))
is.mito <- grepl("^MT-", rownames(sce.covid))

## add mito and spike cols to sce object, which will be present in col-data slot
sce.covid <- addPerCellQC(sce.covid, subsets=list( Mito=is.mito, Spike=is.spike))
## add feature QC data
sce.covid <- addPerFeatureQC(sce.covid)
rowData(sce.covid)

table(as.factor(sce.covid$cluster_id)) ## of cells in smallest group

## for demonstration-purposes, we will first apply a global filter to eliminate all genes without at least a count of 2 in any cell
keep_feature <- 
  rowSums(counts(sce.covid) > 1) > 7 ## smallest group of cells ("Baso")
sce.covid  <- sce.covid[keep_feature, ]

## only keep cells that have at least 20 genes expressed at count of 10 or more
keep_cell <- 
  colSums(counts(sce.covid) > 10) > 20
sce.covid  <- sce.covid[,keep_cell]
min(meta_data$nFeature_RNA)
min(meta_data$nCount_RNA)
## check dimensions
sce.covid ## 10933 19856

## cells with high mitochondria-content are typically associated with low quality reads, 
## so we will compute the outlier-threshold for mito-percentage (usually cells between 5-10% 
## mito are filtered out of the data)

qc.mito <- isOutlier(sce.covid$subsets_Mito_percent,  type="higher") ## percentage threshold val
table(qc.mito) ## 2319 cells are above the threshold
attr(qc.mito, "thresholds") # 6.108 -- this percentage is in line with typical best-practice mito-thresholds

### Repeat code with ERCC (spike-ins)
qc.spike <- isOutlier(sce.covid$subsets_Spike_percent, type="lower")
table(qc.spike) #### ALL FALSE, SO WE WILL LEAVE AS IS

### plot mito content
 plotColData(sce.covid, x="sample_name", y="subsets_Mito_percent",
                          colour_by=
                                as.factor(ifelse(sce.covid$subsets_Mito_percent > attr(qc.mito, "thresholds") , 1, 0)) 
                                          + scale_x_log10()


