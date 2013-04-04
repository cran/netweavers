
generateContrastMatrix <- function(design){
  contrast       <- c();
  index          <- 1
  for (i in 1:(ncol(design)-1)){
    for (j in ((i+1):ncol(design))){
      contrast[index] <- paste(colnames(design)[i],"-",colnames(design)[j])
      index           <- index + 1
    }
  }
  arglist        <- vector("list",length(contrast))
  for (i in 1:length(contrast)){
    arglist[i] <- contrast[i]
  }
  arglist$levels <- design
  do.call("makeContrasts",arglist)
}

DEtest<-function(eset,label,reps=FALSE){
  if (!(label %in% varLabels(eset))){
    cat("\nError: There is no factor '",label,"' that can be used for determining differentially expressed features\n",sep="")
    quit()
  }

  npeaks                    <- nrow(exprs(eset))
  
  ## experimental factor
  fact                      <- as.factor(pData(eset)[,label])
  
  ## linear models in limma
  design                    <- model.matrix(~-1 + fact)
  colnames(design)          <- make.names(levels(fact))

  if (reps){
    ## mixed model
    biolrep <- pData(eset)$sampleID
    corfit  <- duplicateCorrelation(eset,design,ndups=1,block=biolrep)
    cat("correlation between technical replicates:",corfit$consensus,"\n")
    fit     <- lmFit(eset,design=design,block=biolrep,correlation=corfit$consensus)
  }
  else {
    fit     <- lmFit(eset,design=design)
  }
  contrast.matrix           <- generateContrastMatrix(design)
  fit.eb                    <- eBayes(contrasts.fit(fit,contrast.matrix))
  ## F-test
  if (ncol(contrast.matrix)>=2){
    cat(label,"(F-test)\n")
    mt.F.pvalue          <- p.adjust(fit.eb$F.p.value,method="fdr")
    ord                  <- order(fit.eb$F.p.value)
    ## rank of contrast matrix via qr
    resultsTable_Fpvalue <- cbind(fit.eb$genes,F=fit.eb$F,df1=qr(contrast.matrix)$rank,df2=fit.eb$df.residual,
                                  p.Val=fit.eb$F.p.value,adj.P.Val=mt.F.pvalue)[ord,]
    if(any(is.na(resultsTable_Fpvalue$F))){
      cat("\nWarning: nestedF multiple testing scheme not carried out due to missing p-values\n")
    } else {
      results_fit <- decideTests(fit.eb,method="nestedF",adjust.method="fdr")
    }
  }

  resultsTable_limma        <- list()
  
  for (i in 1:ncol(contrast.matrix)){
    cat(label,"(pairwise comparison):",colnames(contrast.matrix)[i],"\n")
    # since df is a constant re-ordering is not necessary
    resultsTable_limma[[i]] <- cbind(topTable(fit.eb, coef=i, adjust.method="fdr", 
                                     number=npeaks, sort.by="p"), df=fit.eb$df.residual)    
  }
  names(resultsTable_limma) <- colnames(contrast.matrix)
  
  if(ncol(contrast.matrix)>=2){
    if(any(is.na(resultsTable_Fpvalue$F))){
      return(list(resultsTable_limma=resultsTable_limma,
                  resultsTable_Fpvalue=resultsTable_Fpvalue))
    } else {
      return(list(resultsTable_limma=resultsTable_limma,
                  resultsTable_Fpvalue=resultsTable_Fpvalue,
                  results_fit=results_fit))
    }
  } 
  else{
    return(list(resultsTable_limma=resultsTable_limma))
  }  
}