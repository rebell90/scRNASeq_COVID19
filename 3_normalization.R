
#############################################################################
################################ NORMALIZATION ##############################
#############################################################################


## recompute library size factors -- keeps normalizated counts in proportion with original counts in respects to lib size
lib.sf.covid <- librarySizeFactors(sce.covid)
summary(lib.sf.covid )
hist(log10(lib.sf.covid ), xlab="Log10[Size factor]", col='purple') ## very slight positive skew

## Lib Size Factors assume  that the DE genes are balanced-
## i.e. any up-regulation on sets of genes will have an equal down-regulation.
## But, in single cell datathis is less likely to be the case due to sparsity in the data.
## This phenomenon is known as composition bias

## normalization methods for SCE data are being explored, given that current methods for bulk-seq do not always translate 
## well to SCE

## normalize by deconvolution:
## pool counts from many cells to increase the size of the counts 
## to accurately calculate size factor estimation. Then, these pool-based size factors
## are deconvolved back to cell-based factors for normalization of each cell's 
## expression file.
## cells are first clustered, then normalization is computed within each cluster, 
## followed by rescaling of size factors to be comparable across all cluster types.
set.seed(100)

clust.covid <- quickCluster(sce.covid) 
table(clust.covid) ## 19
deconv.sf.covid <- calculateSumFactors(sce.covid, cluster=clust.covid)
summary(deconv.sf.covid)

## adjusts for cell-specific deviations and biases
plot(lib.sf.covid, deconv.sf.covid, xlab="Library size factor",
     ylab="Deconvolution size factor", log='xy', pch=16,
     col=as.integer(factor(sce.covid$cluster_id)))
abline(a=0, b=1, col="red")

## add log_counts
set.seed(100) 
sce.covid <- computeSumFactors(sce.covid, cluster=clust.covid)
sce.covid <- logNormCounts(sce.covid)
assayNames(sce.covid)

