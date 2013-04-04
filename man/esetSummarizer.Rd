\name{esetSummarizer}
\alias{esetSummarizer}
\title{
  Aggregate Quantified Peptides Values From an \code{ExpressionSet}
}
\description{
  Summarizes quantified peptide values (e.g. intensity, abundance) in the \code{assayData} slot
  of an \href{http://bioconductor.org/packages/release/bioc/vignettes/Biobase/inst/doc/ExpressionSetIntroduction.pdf}{\code{ExpressionSet}} into a single value per protein.
}
\usage{
esetSummarizer(eset = eset, sumprot = "median")
}
\arguments{
  \item{eset}{
    An \code{ExpressionSet} object containing peptide-level data. 
  }
  \item{sumprot}{
    Function (median or mean) for how the peptide values should be aggregated.
  }
}
\details{
This function calls \code{\link{summarizer}} to roll up peptide-level data to protein-level data. 
The argument \code{eset} should have quantified peptide data in the \code{assayData} slot and 
protein IDs lableled ``Protein'' are required in the \code{featureData} slot. 
Additional variables of annotation are allowed in \code{featureData}, e.g. sequences, but are ignored.

The function returns an \code{ExpressionSet} with the same \code{featureData} that was input but the
\code{assayData} contains the summarized, protein-level data.
 
If the peptide-level data are not already in an \code{ExpressionSet}, see \code{\link{customSummarizer}} 
for summarizing data frames containing quantified values and annotation.
}
\value{
  An \code{ExpressionSet} object containing protein-level data.
}
\seealso{
\code{\link{customSummarizer}}, \code{\link{summarizer}}, \code{\link{pvalueSummarizer}}
}
\examples{
data(vanHoof)
esetSummarizer(ExpSetVH)
}
\keyword{ methods }
