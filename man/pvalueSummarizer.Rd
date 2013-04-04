\name{pvalueSummarizer}
\alias{pvalueSummarizer}
\alias{fisherMethod}
\title{
  Fisher's Method for Multiple Peptides
}
\description{
  \code{pvalueSummarizer} aggregates peptide \emph{p}-values 
  into protein \emph{p}-values using \code{fisherMethod}. 
  \code{fisherMethod} calculates the Fisher's Method \emph{p}-value. 
}
\usage{
pvalueSummarizer(peptide_pvals)
fisherMethod(x)
}
\arguments{
   \item{peptide_pvals}{
    Data frame with at least the following two variables: 
    \code{Protein}, a vector of protein identifiers, and \code{pvalue}, 
    a numeric vector of peptide \emph{p}-values (between 0 and 1).
  }
  \item{x}{
    Numeric vector with values between 0 and 1 representing \emph{p}-values.
  }
}
\details{
The function \code{pvalueSummarizer} calls \code{fisherMethod} for each protein with multiple (peptide)
\emph{p}-values. \code{fisherMethod} returns a single value for a set of \emph{p}-values: 
the chi-squared distribution function on \eqn{-2\sum{\ln(x)}}{-2*sum{log(x)}} with degrees of freedom equal to the length of \code{x}. 
}
\value{
  \code{pvalueSummarizer} returns a data frame of identifiers (\code{Protein}) and the summarized \emph{p}-values (\code{omnibus.P.Val}).
  
  \code{fisherMethod} returns a numeric value between 0 and 1.
}
\references{
  Fisher, R.A. (1925). 
  Statistical Methods for Research Workers. 
  Oliver and Boyd (Edinburgh). 
  ISBN 0-05-002170-2.
}
\seealso{
\code{\link{customSummarizer}}, \code{\link{esetSummarizer}}, \code{\link{summarizer}}
}
\examples{
set.seed(1234)
ids <- sample(letters, 50, replace=TRUE)
pvals <- data.frame(Protein=ids, pvalue=runif(50,0,1))
pvalueSummarizer(pvals)
}
\keyword{ methods }
