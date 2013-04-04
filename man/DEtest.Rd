\name{DEtest}
\alias{DEtest}
\alias{generateContrastMatrix}
\title{
  Differential Expression Testing using Linear Models and Contrasts
}
\description{
  The function \code{DEtest} tests linear models on quantified protein- 
  or peptide-level data. 
  Functions from the package \href{http://bioinf.wehi.edu.au/limma/}{\pkg{limma}} (e.g. \code{contrasts.fit}, \code{eBayes}, \code{topTable})
  are used to calculate \emph{p}-values and other statistics. 
  The function \code{generateContrastMatrix} 
  constructs all combinations of pairwise contrasts from the design matrix constructed in \code{DEtest}.
}
\usage{
DEtest(eset, label, reps = FALSE)
generateContrastMatrix(design)
}
\arguments{
  \item{eset}{
    An \code{ExpressionSet} object containing peptide- or protein-level data. 
  }
  \item{label}{
    Name of the variable subject to testing; 
    should be the name of a column in \code{pData(eset)} 
    (i.e. value from \code{varLabels(eset)}).  
  }
  \item{reps}{
    Logical indicator of whether the dataset contains technical replicates. 
  }
  \item{design}{
    Design matrix with rows corresponding to runs/samples and columns to 
    coefficients to be estimated.
  }
}
\details{
If \code{reps} is TRUE, then the phenoData in \code{eset} must have a column ``sampleID'' and the \pkg{limma} function 
\code{duplicateCorrelation} will estimate the intra-replicate correlation. 

The variable subject to testing must be a factor with at least two levels. If more than 
2 levels are present, then a nestedF multiple testing scheme is carried out.
}
\value{
  For \code{DEtest} a list of: 
  \item{resultsTable_limma}{List of data frames, one for each contrast, that 
  are the output of \code{topTable} applied to the \code{eBayes} models fit.}
  \item{resultsTable_Fpvalue}{Data frame containing (moderated) statistics from an \emph{F}-test. 
  This component is not returned if the variable subject to testing has only 2 factor levels.} 
  
  For \code{generateContrastMatrix} a matrix with columns corresponding to contrasts.
}
\references{
Smyth GK. (2004). Linear models and empirical Bayes methods for assessing differential expression in microarray experiments. 
\emph{Statistical Applications in Genetics and Molecular Biology} 3(1) Article 3.
}
\author{
  Perry D. Moerland <p.d.moerland@amc.uva.nl> and 
  
  Elizabeth A. McClellan <e.mcclellan@erasmusmc.nl>
}
\examples{
data(vanHoof)
de_out <- DEtest(summarizedExpSetVH,label="groups")
str(de_out)
}
\keyword{ methods }
\keyword{ design }
\keyword{ models }