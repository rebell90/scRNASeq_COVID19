

#############################################################################
######################## Dimensionality Reduction ###########################
#############################################################################


# With PCA, its assumed that biological variation is captured in earlier Pcs, 
# and the higher PCs are virtually just technical noise]

set.seed(100) # See below.
sce.covid<- runPCA(sce.covid, subset_row=top_10per) 
reducedDimNames(sce.covid)

dim(reducedDim(sce.covid , "PCA"))

library(BiocSingular)
set.seed(1000)
sce.covid  <- runPCA(sce.covid , subset_row=top_10per, 
                     BSPARAM=RandomParam(), name="IRLBA")
reducedDimNames(sce.covid )

percent.var <- attr(reducedDim(sce.covid), "percentVar")
plot(percent.var, log="y", xlab="PC", ylab="Variance explained (%)")

### we can choose how many PCs , or use heuristics to give us a good "elbow-point"
# Percentage of variance explained is tucked away in the attributes.
chosen.elbow <- PCAtools::findElbowPoint(percent.var)
chosen.elbow  # 5
plot(percent.var, xlab="PC", ylab="Variance explained (%)")
abline(v=chosen.elbow, col="red")

## top 5 percent
set.seed(100) # See below.
top_5per <- getTopHVGs(dec.covid, prop=0.05)
sce.covid <- runPCA(sce.covid, subset_row=top_5per) 
reducedDimNames(sce.covid)

library(BiocSingular)
set.seed(1000)
sce.covid <- runPCA(sce.covid , subset_row=top_5per, 
                     BSPARAM=RandomParam(), name="IRLBA")
reducedDimNames(sce.covid)

percent.var.5p <- attr(reducedDim(sce.covid ), "percentVar")
plot(percent.var.5p, log="y", xlab="PC", ylab="Variance explained (%)")

### we can choose how many PCs , or use heuristics to give us a good "elbow-point"
# Percentage of variance explained is tucked away in the attributes.
chosen.elbow.5p <- PCAtools::findElbowPoint(percent.var.5p)
chosen.elbow.5p  # 5
plot(percent.var.5p, xlab="PC", ylab="Variance explained (%)")
abline(v=chosen.elbow.5p, col="red")

#### top5 percent genes does better

### keep 5 PCS (from elbow threshold) in a new slot
reducedDim(sce.covid, "PCA.elbow") <- reducedDim(sce.covid)[,1:chosen.elbow.5p]
reducedDimNames(sce.covid)


### technical noise
library(scran)
set.seed(111001001)
denoised.covid <- denoisePCA(sce.covid, technical=dec.covid)
ncol(reducedDim(denoised.covid)) #lower bound on # of PCs required to retain all biological variation (5)

## Could also consider subpopulations (that they all differ from eachother)
#  these are usually not known ahead, so a good heuristic approach is to use the number of clusters as a proxy 

## consider a number "d" PCs and only select values that will produce no more than "d + 1" clusters--
## this will ensure that distinct population boundaries between the groups, while also retaining biological signal
## in the later PCs
pcs <- reducedDim(sce.covid, "PCA")
choices <- getClusteredPCs(pcs)
val <- metadata(choices)$chosen

plot(choices$n.pcs, choices$n.clusters,
     xlab="Number of PCs", ylab="Number of clusters")
abline(a=1, b=1, col="red")
abline(v=val, col="grey80", lty=2)
reducedDim(sce.zeisel, "PCA.clust") <- pcs[,1:val]


---------------------------------------
#### visualization
plotReducedDim(sce.covid, dimred="PCA", colour_by="cluster_id")
plotReducedDim(sce.covid, dimred="PCA", ncomponents=4,
               colour_by="cluster_id")

## tsne
set.seed(00101001101)

# runTSNE() stores the t-SNE coordinates in the reducedDims
sce.covid<- runTSNE(sce.covid, dimred="PCA")
plotReducedDim(sce.covid, dimred="TSNE", colour_by="cluster_id")


set.seed(100)
sce.covid <- runTSNE(sce.covid, dimred="PCA", perplexity=5)
out5 <- plotReducedDim(sce.covid, dimred="TSNE",
                       colour_by="cluster_id") + ggtitle("perplexity = 5")

set.seed(100)
sce.covid <- runTSNE(sce.covid, dimred="PCA", perplexity=20)
out20 <- plotReducedDim(sce.covid, dimred="TSNE",
                        colour_by="cluster_id") + ggtitle("perplexity = 20")

gridExtra::grid.arrange(out5, out20,  ncol=2)
#### UMAP
set.seed(1100101001)
sce.covid <- runUMAP(sce.covid, dimred="PCA")
plotReducedDim(sce.covid, dimred="UMAP", colour_by="cluster_id")





