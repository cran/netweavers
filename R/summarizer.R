summarizer <- function(exprs,fdata,sumprot="median"){
  if(!sumprot%in%c("mean","median")){
     cat("\nError: The argument sumprot must be 'mean' or 'median'.\n",sep="")
     quit()
  }
  n                <- ncol(exprs)
  sumpp            <- apply(exprs,2,function(x){aggregate(x~fdata$Protein,FUN=sumprot)})
  newsum           <- do.call('rbind',lapply(1:length(sumpp),function(x){cbind(sumpp[[x]], q=x)}))
  sumup            <- reshape(newsum, v.names='x', idvar='fdata$Protein', timevar='q',direction='wide')
  colnames(sumup)  <- c("Protein", colnames(exprs))
  row.names(sumup) <- 1:nrow(sumup)
  sumup
}

customSummarizer <- function(pepquant,samplecols,peptidecol,proteincol,sumprot="median",phenodata){
  if(!sumprot%in%c("mean","median")){
     cat("\nError: The argument sumprot must be 'mean' or 'median'.\n",sep="")
     quit()
  }
  if(!all(names(pepquant[,samplecols])==phenodata$sampleIDs)){
     cat("\nError: The names of pepquant[,samplecols] must be identical to and in the same order as phenodata[,'sampleIDs'].\n",sep="")
     quit()
  }
  options(warn=-1)
  ## check parameters
  if(missing(proteincol)) stop("protein IDs are missing") 
  if(missing(samplecols)) stop("sample columns are missing") 
  ## get peptide intensities and feature data
  exprs_x        <- pepquant[,samplecols] ## make sure no 0s in dataset! 
  fdata_x        <- pepquant[,c(proteincol,peptidecol)]
  names(fdata_x) <- c("Protein","Sequence")
  
  sumup          <- summarizer(exprs_x,fdata_x,sumprot)

  ## new ExpressionSet where proteins are features instead of peptides/peaks
  vardf          <- data.frame(names=names(phenodata),labelDescription=NA)
  new("ExpressionSet", 
      exprs = as.matrix(sumup[,-1]), 
      phenoData = new("AnnotatedDataFrame",data=phenodata,varMetadata=vardf), 
      featureData = new("AnnotatedDataFrame",data=data.frame(Protein=sumup[,'Protein'])))
}

esetSummarizer <- function(eset=eset,sumprot="median"){
  if(!sumprot%in%c("mean","median")){
     cat("\nError: The argument sumprot must be 'mean' or 'median'.\n",sep="")
     quit()
  }
  ## get expression data.frame
  fdata_x  <- fData(eset)
  exprs_x  <- exprs(eset)

  sumup    <- summarizer(exprs_x,fdata_x,sumprot)

  ## new ExpressionSet where proteins are features instead of peptides/peaks
  new("ExpressionSet", 
      exprs = as.matrix(sumup[,-1]), 
      phenoData = phenoData(eset), 
      featureData = new("AnnotatedDataFrame",data=data.frame(Protein=sumup[,'Protein'])))
}

fisherMethod <- function(x){
  y <- (-2) * sum(log(x))
  pchisq(y, df=length(x),lower.tail=F)
}

pvalueSummarizer <- function(peptide_pvals){ 
  fishtest           <- aggregate(pvalue~Protein,data=peptide_pvals,FUN=fisherMethod)
  names(fishtest)[2] <- "omnibus.P.Val"
  if(class(fishtest$Protein)=="factor"){
    fishtest$Protein   <- levels(fishtest$Protein)[fishtest$Protein]
  }
  fishtest
}