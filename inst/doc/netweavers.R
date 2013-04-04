### R code from vignette source 'netweavers.Rnw'

###################################################
### code chunk number 1: netweavers.Rnw:38-39
###################################################
library(netweavers)


###################################################
### code chunk number 2: netweavers.Rnw:47-48
###################################################
data(vanHoof)


###################################################
### code chunk number 3: netweavers.Rnw:61-62
###################################################
ExpSetVH


###################################################
### code chunk number 4: netweavers.Rnw:67-68
###################################################
phenoData(ExpSetVH)@data


###################################################
### code chunk number 5: netweavers.Rnw:73-74
###################################################
outDE <- DEtest(ExpSetVH, label="groups")


###################################################
### code chunk number 6: netweavers.Rnw:87-88
###################################################
str(outDE)


###################################################
### code chunk number 7: netweavers.Rnw:92-93
###################################################
outDE_30to60 <- outDE$resultsTable_limma[["min30 - min60"]]


###################################################
### code chunk number 8: netweavers.Rnw:106-108
###################################################
head(dataVH)
phenoDataVH


###################################################
### code chunk number 9: netweavers.Rnw:112-116
###################################################
esetSum1 <- customSummarizer(pepquant=dataVH, samplecols=3:8, 
                            peptidecol=1, proteincol=2, 
                            sumprot="median", phenodata=phenoDataVH)
esetSum1


###################################################
### code chunk number 10: netweavers.Rnw:133-134
###################################################
esetSum2 <- esetSummarizer(ExpSetVH, sumprot="median")


###################################################
### code chunk number 11: netweavers.Rnw:152-156
###################################################
pvals3060 <- outDE_30to60[,c("Protein","P.Value")]
names(pvals3060)[2] <- "pvalue"
pvalsSum <- pvalueSummarizer(pvals3060)
head(pvalsSum)


###################################################
### code chunk number 12: netweavers.Rnw:171-172
###################################################
filtnet <- filterNetwork(networkVH, pvalsSum$Protein)


###################################################
### code chunk number 13: netweavers.Rnw:177-182
###################################################
## note: this will take a couple of minutes
## see clust_2_proteinVH for identical output
outDC <- findDenseClusters(filtnet, pvalsSum$Protein, min_clus_size=4, 
                           steps=10)
denseClusters <- outDC$clust_2_protein


###################################################
### code chunk number 14: netweavers.Rnw:188-197
###################################################
## change name of p-values for use in measureCalc
names(pvalsSum)[2] <- "pvalue"

## calculate protein weights
proteinWts <- weightCalc(filtnet)

## calculate protein measures
proteinMsr  <- measureCalc(pvalsSum, proteinWts, weightamt=3)
head(proteinMsr)


###################################################
### code chunk number 15: netweavers.Rnw:202-203
###################################################
cscores <- scoreClusters(denseClusters, proteinMsr, cscoremethod="mean")


###################################################
### code chunk number 16: netweavers.Rnw:208-212
###################################################
set.seed(1234)
cpvals <- permTest(B=1000, prot_data=proteinMsr, 
                   clust_2_protein=denseClusters,
                   cluster_scores=cscores, cscoremethod="mean")


###################################################
### code chunk number 17: netweavers.Rnw:216-221
###################################################
cdata <- data.frame(Cluster_ID=denseClusters$clusID,
                    Cluster_size=denseClusters$clusSIZE,
                    Cluster_score=cscores,
                    Cluster_pvalue=cpvals,
                    Cluster_symbol=denseClusters[,-c(1,2)])


###################################################
### code chunk number 18: netweavers.Rnw:228-231
###################################################
set.seed(1234)
outNW <- runNetweavers(pvalsSum, filtnet, weightamt=3)
str(outNW,list.len=5)


###################################################
### code chunk number 19: netweavers.Rnw:238-244
###################################################
proteininfo <- outNW$proteininfo
clusterinfo <- outNW$clusterinfo
proteininfo$cluster <- NA
for(i in 1:nrow(clusterinfo)){
  proteininfo$cluster[proteininfo[,1]%in%clusterinfo[i,-c(1:4)]] <- i
}


###################################################
### code chunk number 20: netweavers.Rnw:248-250
###################################################
## write.table(proteininfo,"nodeAttributes.txt",
##             sep="\t",row.names=F,quote=F,na="-9")


###################################################
### code chunk number 21: netweavers.Rnw:253-254
###################################################
## write.table(filtnet,"network.txt",sep="\t",row.names=F,quote=F)


###################################################
### code chunk number 22: netweavers.Rnw:269-276
###################################################
dc_member <- outDC[['denseClus']]$membership
dc_graph  <- outDC$subGraph
cols      <- rainbow(length(unique(dc_member)))

## make the plot
## plot(dc_graph, vertex.color=cols, vertex.size=3, vertex.label=NA, 
##      edge.color="grey")


