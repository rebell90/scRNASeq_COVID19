############# marker gene detection

markers.covid <- findMarkers(covid.sce)
markers.covid

chosen <- "5"
interesting <- markers.covid[[chosen]]
colnames(interesting)

best.set <- interesting[interesting$Top <= 5,]
logFCs <- getMarkerEffects(best.set)

library(pheatmap)
pheatmap(logFCs, breaks=seq(-5, 5, length.out=101))

### DE genes in cluster of interest compared to all other clsuter

# Set direction='up' to only consider upregulated genes as potential markers.
markers.covid.up3 <- findMarkers(covid.sce, pval.type="all", direction="up")
interesting.up3 <- markers.covid.up3[[chosen]]
interesting.up3[1:10,1:3]

markers.covid.up4 <- findMarkers(covid.sce, pval.type="some", direction="up")
interesting.up4 <- markers.covid.up4[[chosen]]
interesting.up4[1:10,1:3]


markers.covid.up <- findMarkers(covid.sce, direction="up")
interesting.up <- markers.covid.up[[chosen]]
interesting.up[1:10,1:4]

markers.covid.up2 <- findMarkers(covid.sce, direction="up", lfc=1)
interesting.up2 <- markers.covid.up2[[chosen]]
interesting.up2[1:10,1:4]

best.set <- interesting.up2[interesting.up2$Top <= 5,]
logFCs <- getMarkerEffects(best.set)

library(pheatmap)
pheatmap(logFCs, breaks=seq(-5, 5, length.out=101))


## wilcox rank
markers.covid.wmw <- findMarkers(covid.sce, test="wilcox", direction="up")
names(markers.covid.wmw)

interesting.wmw <- markers.covid.wmw[[chosen]]
interesting.wmw[1:10,1:4]

best.set <- interesting.wmw[interesting.wmw$Top <= 5,]
AUCs <- getMarkerEffects(best.set, prefix="AUC")

library(pheatmap)
pheatmap(AUCs, breaks=seq(0, 1, length.out=21),
         color=viridis::viridis(21))

------------
  marker.covid.t <- findMarkers(covid.sce, groups=covid.sce$`cell type`, 
                                 direction="up", restrict=c("Alpha", "Gamma/PP"))
marker.covid.w <- findMarkers(covid.sce, groups=covid.sce$`cell type`, 
                               direction="up", restrict=c("Alpha", "Gamma/PP"), test.type="wilcox")

# Upregulated in alpha:
marker.alpha.t <- marker.covid.t$Alpha
marker.alpha.w <- marker.covid.w$Alpha
chosen.alpha.t <- rownames(marker.alpha.t)[1:20]
chosen.alpha.w <- rownames(marker.alpha.w)[1:20]



### binomial test
markers.covid.binom <- findMarkers(covid.sce, test="binom", direction="up")
names(markers.covid.binom)

interesting.binom <- markers.covid.binom[[chosen]]
colnames(interesting.binom)

library(scater)
top.genes <- head(rownames(interesting.binom))
plotExpression(covid.sce, x="label", features=top.genes)

plotExpression(covid.sce, x="cell type", features="ENSG00000115263")
u.alpha.t <- setdiff(chosen.alpha.t, chosen.alpha.w)
u.alpha.w <- setdiff(chosen.alpha.w, chosen.alpha.t)

# Upregulated in gamma:
marker.gamma.t <- marker.covid.t$`Gamma/PP`
marker.gamma.w <- marker.covid.w$`Gamma/PP`
chosen.gamma.t <- rownames(marker.gamma.t)[1:20]
chosen.gamma.w <- rownames(marker.gamma.w)[1:20]
u.gamma.t <- setdiff(chosen.gamma.t, chosen.gamma.w)
u.gamma.w <- setdiff(chosen.gamma.w, chosen.gamma.t)

# Examining all uniquely detected markers in each direction.
library(scater)
subset <- covid.sce[,covid.sce$`cell type` %in% c("Alpha", "Gamma/PP")]
gridExtra::grid.arrange(
  plotExpression(subset, x="cell type", features=u.alpha.t, ncol=2) +
    ggtitle("Upregulated in alpha, t-test-only"),
  plotExpression(subset, x="cell type", features=u.alpha.w, ncol=2) +
    ggtitle("Upregulated in alpha, WMW-test-only"),
  plotExpression(subset, x="cell type", features=u.gamma.t, ncol=2) +
    ggtitle("Upregulated in gamma, t-test-only"),
  plotExpression(subset, x="cell type", features=u.gamma.w, ncol=2) +
    ggtitle("Upregulated in gamma, WMW-test-only"),
  ncol=2
)

