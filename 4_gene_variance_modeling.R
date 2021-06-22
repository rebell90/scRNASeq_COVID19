#############################################################################
################################ Variance Modeling ##############################
#############################################################################

#### calc within gene variance--- traditional methods to calc mean-variance relationship
## have underlying assumptions on heteroscedasciity , which SCE data rarely has, thus thed values are driven more by the 
## abundance---  using modelGeneVar() function will help account for this 
dec.covid <- modelGeneVar(sce.covid)

plot(dec.covid$mean, dec.covid$var, xlab="Mean of log-expression",
     ylab="Variance of log-expression")


# Ordering by most interesting genes
dec.covid[order(dec.covid$bio, decreasing=TRUE),] 

## most sig genes
hvg.covid.var <- getTopHVGs(dec.covid, n=1000)
str(hvg.covid.var)
# sig threshold
hvg.covid.var.2 <- getTopHVGs(dec.covid, fdr.threshold=0.05)
length(hvg.covid.var.2) #140
##above trend
hvg.covid.var.3 <- getTopHVGs(dec.covid, var.threshold=0)
length(hvg.covid.var.3) #7434

## top 10 percent of genes with highest biological components
top_10per <- getTopHVGs(dec.covid, prop=0.1)
str(top_10per)
length(top_10per)

## coefficient of variation
dec.cv2.covid <- modelGeneCV2(sce.covid)
hvg.covid.cv2 <- getTopHVGs(dec.cv2.covid, var.field="ratio", n=1000)
str(hvg.covid.cv2)
hvg.covid.cv2.3 <- getTopHVGs(dec.cv2.covid, var.field="ratio", var.threshold=1)
length(hvg.covid.cv2.3) # 8246

###### quanitify technical noise ( in absence of spike in data)
set.seed(0010101)
dec.pois.covid<- scran::modelGeneVarByPoisson(sce.covid)
dec.pois.covid <- dec.pois.covid[order(dec.pois.covid$bio, decreasing=TRUE),]
head(dec.pois.covid)

## mean-variance of log-expression
plot(dec.pois.covid$mean, dec.pois.covid$total, pch=16, xlab="Mean of log-expression",
     ylab="Variance of log-expression")
## this method usually yields large biological components for 
## highly expressed genes, even those "house-keeping" genes which are normally assumed to be un-interesting, 
## but it is still important to look/analyze this to see if it could lead to any unqiue/interesting insights

# Performing PCA only on the chosen HVGs.
library(scater)
sce.covid <- runPCA(sce.covid, subset_row=top_10per)
reducedDimNames(sce.covid )


####### store objects in SCE
rowSubset(sce.covid) <- "top_10per"  # stored as default ( calculated above)
rowSubset(sce.covid, "top_5per") <- getTopHVGs(dec.covid, prop=0.05)
rowSubset(sce.covid, "top_20per") <- getTopHVGs(dec.covid, prop=0.2)
rowSubset(sce.covid, "top_30per") <- getTopHVGs(dec.covid, prop=0.3)
rowSubset(sce.covid, "cv2_top_1000") <- hvg.covid.cv2 
rowSubset(sce.covid, "cv2_top_var_thresh") <- hvg.covid.cv2.3 

colnames(rowData(sce.covid))
str(top_10per)
