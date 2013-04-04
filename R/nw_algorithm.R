
fillup <- function(mylist){
  newcol  <- max(sapply(mylist,length))
  newlist <- lapply(mylist,function(x){c(x,rep(NA,newcol))[1:newcol]})
  do.call(rbind,newlist)
}

filterNetwork <- function(network,identifiers){
  filt_net1 <- network[network[,1]%in%identifiers | network[,2]%in%identifiers,]
  unique(filt_net1[!duplicated(t(apply(filt_net1,1,sort))),])
}

findDenseClusters <- function(filt_net,identifiers,min_clus_size=4,steps=10){
  subGraph <- simplify(graph.edgelist(unique(as.matrix(filt_net)),directed=F))
  nodeProtein <- get.vertex.attribute(subGraph,"name")
  denseClus <- walktrap.community(subGraph,steps=steps)
  dc1 <- denseClus$member
  names(dc1) <- nodeProtein
  dc2 <- list()
  for(i in 1:max(dc1)){dc2[[i]] <- names(dc1[dc1==i])}
  dc3             <- lapply(dc2,function(x){
                                  if(length(x)<min_clus_size) x<-NA else x<-x
                                  if(all(!x%in%identifiers)) x<-NA else x<-x
                                  })
  dc4             <- dc3[!is.na(dc3)]
  clust_2_protein <- data.frame(clusID=c(paste("c",1:length(dc4),sep="")),
                                clusSIZE=sapply(dc4,length),fillup(dc4))
  ## do not want factors - convert to character
  for(i in 3:ncol(clust_2_protein)){
      clust_2_protein[,i] <- levels(clust_2_protein[,i])[clust_2_protein[,i]]
    }
  list(clust_2_protein=clust_2_protein,subGraph=subGraph,denseClus=denseClus)
}

weightCalc <- function(filt_net){  
  all_gene_int                  <- merge(table(filt_net[,2]),table(filt_net[,1]),all=TRUE)
  neighbors                     <- aggregate(all_gene_int[,2],by=list(all_gene_int[,1]),sum)
  names(neighbors)              <- c("Protein","n")
  #neighbors[neighbors$n==0,'n'] <- NA ## 0 shouldn't be possible, but this occurs with extraneous factor levels
  data.frame(Protein=neighbors$Protein,weight=1/neighbors$n)
}

measureCalc <- function(protein_pvals,protein_wts,weightamt = 10){
  prot_data                        <- merge(protein_pvals,protein_wts,by='Protein',all=TRUE)
  prot_data$measure                <- -log((prot_data$pvalue)*(prot_data$weight)^(1/weightamt))
  net_no_pval                      <- is.na(prot_data$measure) & !is.na(prot_data$weight)
  prot_data[net_no_pval,'measure'] <- 0
  prot_data
}

scoreClusters <- function(clust_2_protein,prot_meas,cscoremethod="mean"){
  if(!cscoremethod%in%c("mean","median")){
     cat("\nError: The argument cscoremethod must be 'mean' or 'median'.\n",sep="")
     quit()
  }
  clust_score_total        <- apply(clust_2_protein[,-c(1,2)],1,
                                    function(x){do.call(cscoremethod,list(prot_meas[prot_meas$Protein%in%x[!is.na(x)],'measure']))})
  names(clust_score_total) <- clust_2_protein[,1]
  clust_score_total
}

permScores <- function(perm_prot,prot_data,clust_2_protein,cscoremethod="mean"){
  if(!cscoremethod%in%c("mean","median")){
     cat("\nError: The argument cscoremethod must be 'mean' or 'median'.\n",sep="")
     quit()
  }
  new_prot_data         <- prot_data
  new_prot_data$Protein <- perm_prot
  scoreClusters(clust_2_protein,new_prot_data,cscoremethod)
}

permTest <- function(B,prot_data,clust_2_protein,cluster_scores,cscoremethod="mean"){ 
  if(!cscoremethod%in%c("mean","median")){
     cat("\nError: The argument cscoremethod must be 'mean' or 'median'.\n",sep="")
     quit()
  }
  B                <- as.numeric(B)
  prot_data        <- prot_data[!is.na(prot_data$weight),] ## only permute proteins in network
  perms            <- data.frame(og=prot_data[,1]) 
  for(b in 1:B){
    perms <- data.frame(perms,b=sample(perms[,1]))
  }
  all_cscore_perms <- apply(perms[,-1],2,permScores,prot_data,clust_2_protein,cscoremethod)
  new_forprops     <- cbind(cluster_scores,all_cscore_perms)
  apply(new_forprops,1,function(x){(length(x[x>=x[1]])-1)/B})
}
   


runNetweavers<-function(protein_pvals,network,clust_2_protein=NA,min_clus_size=4,
  steps=10,weightamt=10,cscoremethod="mean",permtest=TRUE,B=1000){
  
  if(any(!(c("Protein","pvalue") %in% names(protein_pvals)))){
     cat("\nError: The argument protein_pvals must contain variables 'Protein' and 'pvalue'.\n",sep="")
     quit()
  }
  if(!cscoremethod%in%c("mean","median")){
     cat("\nError: The argument cscoremethod must be 'mean' or 'median'.\n",sep="")
     quit()
  }
  
  cat("filtering network\n")
  filt_net          <- filterNetwork(network,protein_pvals$Protein)

  if(missing(clust_2_protein)){
    cat("finding dense clusters\n")
    clust_2_protein <- findDenseClusters(filt_net,protein_pvals$Protein,
	                                     min_clus_size,steps)$clust_2_protein
  }

  cat("calculating proteins weights\n")
  protein_wts        <- weightCalc(filt_net)

  cat("calculating protein measures\n")
  prot_data         <- measureCalc(protein_pvals,protein_wts,weightamt)
  
  cat("calculating cluster scores\n")
  cluster_scores    <- scoreClusters(clust_2_protein,prot_data[,c('Protein','measure')],cscoremethod)

  if(permtest==TRUE){
    cat("calculating permutation pvalues\n")
    cluster_pvals   <- permTest(B,prot_data,clust_2_protein,cluster_scores,cscoremethod)
  } else{ 
      cluster_pvals <- NA 
      }
  
  cat("producing results tables\n")
  clusterinfo       <- data.frame(Cluster_ID=clust_2_protein$clusID,
                                  Cluster_size=clust_2_protein$clusSIZE,
                                  Cluster_score=cluster_scores,
                                  Cluster_pvalue=cluster_pvals,
                                  Cluster_symbol=clust_2_protein[,-c(1,2)])

  list(proteininfo=prot_data,filt_net=filt_net,clusterinfo=clusterinfo)
}























