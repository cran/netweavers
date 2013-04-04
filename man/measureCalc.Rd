\name{measureCalc}
\alias{measureCalc}
\title{
  Calculate Protein Measures for NetWeAvers Algorithm
}
\description{
  Calculates a weighted measure for each protein in an interaction network. 
}
\usage{
measureCalc(protein_pvals, protein_wts, weightamt = 10)
}
\arguments{
  \item{protein_pvals}{
    Data frame on two variables: 
    \code{Protein}, a factor or vector of protein identifiers, and \code{pvalue}, 
    a numeric vector of protein \emph{p}-values (between 0 and 1).
  }
  \item{protein_wts}{
    Data frame on two variables: 
    \code{Protein}, a factor or vector of protein identifiers, and \code{weight}, 
    a numeric vector of protein weights (between 0 and 1).
  }
  \item{weightamt}{
    Positive numeric value (see Details). 
  }
}
\details{
This function is meant for use in the NetWeAvers algorithm, where the weight 
for a given protein is the inverse of the number of proteins with which it interacts. 

The argument \code{weightamt} should be a number indicating how much to downweight the protein weights as compared to the protein \emph{p}-values. 
The larger this value, the larger the impact of \emph{p}-values on the protein measure as compared to the weights.
When \code{weightamt} is 1, the \emph{p}-values are equally as important as the weights; the closer it is to 0, the more important the weights. 
The default \code{weightamt} is 10, which makes the \emph{p}-values more important than the weights in the algorithm, as does any number larger than 1.  

The equation for calculating the protein measure is \eqn{-\ln(pvalue*weight^{(1/weightamt)})}{-ln(pvalue*weight^(1/weightamt))}. 
The two data frame arguments are merged on the "Protein" variable in each data frame, so at least a subset of the identifiers should match.  
A measure of 0 is given to all proteins with weights but no \emph{p}-values.
}
\value{
  Data frame on the following 4 variables:
  \describe{
    \item{\code{Protein}}{Factor or vector of protein identifiers; the result of merging \code{protein_pvals} and \code{protein_wts}.}
    \item{\code{pvalue }}{Numeric vector of supplied \emph{p}-values; \code{NA} where weight present but \emph{p}-value not.}
    \item{\code{weight }}{Numeric vector of supplied weights; \code{NA} where \emph{p}-value present but weight not.}
    \item{\code{measure}}{Numeric vector of values calculated as described above (see Details).}
  }
}
\seealso{
\code{\link{weightCalc}}, \code{\link{runNetweavers}}
}
\examples{
set.seed(1234)
ids <- letters[1:10]
pvalues <- data.frame(Protein = ids[-2], pvalue = runif(9,0,1))
protein_wts <- data.frame(Protein = ids[-4], weight = abs(runif(9)))
measureCalc(pvalues, protein_wts)
}
\keyword{ methods }
