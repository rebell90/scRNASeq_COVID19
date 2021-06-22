

#############################################################################
############################### CLUSTERING ##################################
#############################################################################

library(scran)
g <- buildSNNGraph(sce.covid, k=10, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
table(clust)


library(bluster)
clust2 <- clusterRows(reducedDim(sce.covid, "PCA"), NNGraphParam())
table(clust2) # same as above.

library(scater)
colLabels(sce.covid) <- factor(clust)
plotReducedDim(sce.covid, "TSNE", colour_by="label")

# More resolved.
g.5 <- buildSNNGraph(sce.covid, k=5, use.dimred = 'PCA')
clust.5 <- igraph::cluster_walktrap(g.5)$membership
table(clust.5)

# Less resolved.
g.50 <- buildSNNGraph(sce.covid, k=50, use.dimred = 'PCA')
clust.50 <- igraph::cluster_walktrap(g.50)$membership
table(clust.50)

set.seed(11000)
reducedDim(sce.covid, "force") <- igraph::layout_with_fr(g)
plotReducedDim(sce.covid, colour_by="label", dimred="force")

## weighting parameters
## number weights edges based on the # of nearest neighbors that are shared between two cells
## jaccard uses the jaccard index to weight
## none means no weighting

g.num <- buildSNNGraph(sce.covid, use.dimred="PCA", type="number")
g.jaccard <- buildSNNGraph(sce.covid, use.dimred="PCA", type="jaccard")
g.none <- buildKNNGraph(sce.covid, use.dimred="PCA")

## other algorithms for constructing graph object
clust.louvain <- igraph::cluster_louvain(g)$membership
clust.infomap <- igraph::cluster_infomap(g)$membership
clust.fast <- igraph::cluster_fast_greedy(g)$membership
clust.labprop <- igraph::cluster_label_prop(g)$membership
clust.eigen <- igraph::cluster_leading_eigen(g)$membership


library(pheatmap)

# Using a large pseudo-count for a smoother color transition
# between 0 and 1 cell in each 'tab'.
tab <- table(paste("Infomap", clust.infomap), 
             paste("Walktrap", clust))
ivw <- pheatmap(log10(tab+10), main="Infomap vs Walktrap",
                color=viridis::viridis(100), silent=TRUE)

tab <- table(paste("Fast", clust.fast), 
             paste("Walktrap", clust))
fvw <- pheatmap(log10(tab+10), main="Fast-greedy vs Walktrap",
                color=viridis::viridis(100), silent=TRUE)

gridExtra::grid.arrange(ivw[[4]], fvw[[4]])

# infomap > walktrap > fast-greedy (# of clusters)


## hierarchal (combining clusters on similarity)
community.walktrap <- igraph::cluster_walktrap(g)
table(igraph::cut_at(community.walktrap, n=5)) ## very uneven
## try n = 20
table(igraph::cut_at(community.walktrap, n=20)) ## a little better

#for non-heirachal , mergeCommunities() functioin will greedily merge pairs of clusters 
# until a specified number of clusters is achieved

community.louvain <- igraph::cluster_louvain(g)
table(community.louvain$membership)
merged <- mergeCommunities(g, community.louvain$membership, number=10)
table(merged) ### much more even ! 


## analyzing seperation between clusters
# using the ratio of the observed to expected sum of weights between each pair of clusters.
## this will give a better result than the "difference" , because it will be less dependent on population per cluster
ratio <- pairwiseModularity(g, clust, as.ratio=TRUE)
dim(ratio)

#each row/column corresponds to a cluster and entry is the ratios of observed to total weight of edges between cells in the respective clusters. 
## well seperated clusters should have higher weight on diagonal because most edges occur between cells in the same cluster
library(pheatmap)
pheatmap(log2(ratio+1), cluster_rows=FALSE, cluster_cols=FALSE,
         color=colorRampPalette(c("white", "blue"))(100))


## cluster graph to explore relationships between clusters (can be used for exploration)
cluster.gr <- igraph::graph_from_adjacency_matrix(log2(ratio+1), 
                                                  mode="upper", weighted=TRUE, diag=FALSE)

# Increasing the weight to increase the visibility of the lines.
set.seed(11001010)
plot(cluster.gr, edge.width=igraph::E(cluster.gr)$weight*5,
     layout=igraph::layout_with_lgl)

### K-means and heirachal too computationially expensive
######### Cluster Diagnostics #########
## cluster sepeation
# Performing the calculations on the PC coordinates, like before.

#average of the distances with the distance to the average (i.e., centroid) of each cluster, 
#with some tweaks to account for the distance due to the within-cluster variance.
sil.approx <- approxSilhouette(reducedDim(sce.covid, "PCA"), clusters=clust)

sil.data <- as.data.frame(sil.approx)
sil.data$closest <- factor(ifelse(sil.data$width > 0, clust, sil.data$other))
sil.data$cluster <- factor(clust)

ggplot(sil.data, aes(x=cluster, y=width, colour=closest)) + 
  ggbeeswarm::geom_quasirandom(method="smiley")

## cluster purity -- are how much do the cells intertwine/mingle together in the expression space?\
## well defined clusters should have little to no instances of this
pure.covid <- neighborPurity(reducedDim(sce.covid, "PCA"), clusters=clust)

pure.data <- as.data.frame(pure.covid)
pure.data$maximum <- factor(pure.data$maximum)
pure.data$cluster <- factor(clust)

ggplot(pure.data, aes(x=cluster, y=purity, colour=maximum)) +
  ggbeeswarm::geom_quasirandom(method="smiley")

# purity does not consider intra-cluster variance
##even if clusters are not pure, this is not necessarily a bad thing, it all depends on the situation/data


### different clustering algorithms often lead to different results, but this step is
### merely exploratory and can give us cues to what may need further analysis , so comparing the 
### various clustering algorithms can provide insights -- how are the clusters related in each? are some cells
### being clustered in one group for some algorithm, but in another group with a different algorithm?

#capture the redistribution of cells from one clustering to another
install.packages("clustree")
library(clustree)
combined <- cbind(k.50=clust.50, k.10=clust, k.5=clust.5)
clustree(combined, prefix="k.", edge_arrow=FALSE)

## ARI ( the closer to 1, the more defined the clusters)
pairwiseRand(clust, clust.5, mode="index") ##.88 



breakdown <- pairwiseRand(ref=clust, alt=clust.5, mode="ratio") 
pheatmap(breakdown, color=viridis::magma(100), 
         cluster_rows=FALSE, cluster_cols=FALSE)


## Cluster stability
myClusterFUN <- function(x) {
  g <- bluster::makeSNNGraph(x, type="jaccard")
  igraph::cluster_louvain(g)$membership
}

pcs <- reducedDim(sce.covid, "PCA")
originals <- myClusterFUN(pcs)
table(originals) # inspecting the cluster sizes.

set.seed(0010010100)
ratios <- bootstrapStability(pcs, FUN=myClusterFUN, clusters=originals)
dim(ratios)

pheatmap(ratios, cluster_row=FALSE, cluster_col=FALSE,
         color=viridis::magma(100), breaks=seq(-1, 1, length.out=101))

g.full <- buildSNNGraph(sce.covid, use.dimred = 'PCA')
clust.full <- igraph::cluster_walktrap(g.full)$membership
plotExpression(sce.covid, features=c("CD3E", "CCR7", "CD69", "CD44"),
               x=I(factor(clust.full)), colour_by=I(factor(clust.full)))

# Repeating modelling and PCA on the subset.
memory <- 10L
sce.memory <- sce.covid[,clust.full==memory]
dec.memory <- modelGeneVar(sce.memory)
sce.memory <- denoisePCA(sce.memory, technical=dec.memory,
                         subset.row=getTopHVGs(dec.memory, n=5000))

g.memory <- buildSNNGraph(sce.memory, use.dimred="PCA")
clust.memory <- igraph::cluster_walktrap(g.memory)$membership
plotExpression(sce.memory, features=c("CD8A", "CD4"),
               x=I(factor(clust.memory)))

set.seed(1000010)
subcluster.out <- quickSubCluster(sce.covid, groups=clust.full,
                                  prepFUN=function(x) { # Preparing the subsetted SCE for clustering.
                                    dec <- modelGeneVar(x)
                                    input <- denoisePCA(x, technical=dec,
                                                        subset.row=getTopHVGs(dec, prop=0.1),
                                                        BSPARAM=BiocSingular::IrlbaParam())
                                  },
                                  clusterFUN=function(x) { # Performing the subclustering in the subset.
                                    g <- buildSNNGraph(x, use.dimred="PCA", k=20)
                                    igraph::cluster_walktrap(g)$membership
                                  }
)

# One SingleCellExperiment object per parent cluster:
names(subcluster.out)

# Looking at the subclustering for one example:
table(subcluster.out[[1]]$subcluster)


