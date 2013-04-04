\name{customSummarizer}
\alias{customSummarizer}
\title{
Aggregate Quantified Peptides Values From a Data Frame
}
\description{
  Summarizes quantified peptide values (e.g. intensity, abundance) into a
  single value per protein and returns results in an \code{ExpressionSet}. 
}
\usage{
customSummarizer(pepquant, samplecols, peptidecol, proteincol, 
                 sumprot = "median", phenodata)
}
\arguments{
  \item{pepquant}{
    Data frame of quantified peptide values.
  }
  \item{samplecols}{
    Vector of integers indicating which columns in \code{pepquant} should be summarized.
  }
  \item{peptidecol}{
    Integer; column index in \code{pepquant} for peptide identifiers.
  }
  \item{proteincol}{
    Integer; column index in \code{pepquant} for protein identifiers.
  }
  \item{sumprot}{
    Function (median or mean) for how the peptide values should be aggregated. 
  }
  \item{phenodata}{
    Data frame of phenotypic data (e.g. covariates), 
    including a variable \code{sampleIDs} that is a vector identical to the names of \code{pepquant[,samplecols]}. 
  }
}
\details{
This function calls \code{\link{summarizer}} to roll up peptide-level data to protein-level data. 
Note that in \code{pepquant} the samples are stored column-wise and in \code{phenodata} the samples are stored row-wise. 
Additional variables of annotation, e.g. sequences, are allowed in \code{pepquant} but are ignored.

The function returns an \href{http://bioconductor.org/packages/release/bioc/vignettes/Biobase/inst/doc/ExpressionSetIntroduction.pdf}{\code{ExpressionSet}}
with the summarized, protein-level data in the \code{assayData} slot.
The slot \code{phenoData} contains the \code{phenodata} and the slot \code{featureData}
contains the protein IDs from \code{proteincol}.
 
If the peptide-level data are already in an \code{ExpressionSet}, see \code{\link{esetSummarizer}}.
}
\value{
  An \code{ExpressionSet} containing protein-level data.
}
\seealso{
\code{\link{esetSummarizer}}, \code{\link{summarizer}}, \code{\link{pvalueSummarizer}}
}
\examples{
data(vanHoof)
customSummarizer(dataVH, samplecols=3:8, peptidecol=1, 
                 proteincol=2, phenodata=phenoDataVH)
}
\keyword{ methods }
