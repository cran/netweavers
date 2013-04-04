\name{summarizer}
\alias{summarizer}
\title{
  Aggregate Quantified Peptides Values
}
\description{
  Summarizes quantified peptide values (e.g. intensity, abundance) into a
  single value per protein. 
}
\usage{
summarizer(exprs, fdata, sumprot = "median")
}
\arguments{
  \item{exprs}{
    Data frame of quantified peptide values.
  }
  \item{fdata}{
    Matrix or data frame containing protein IDs and (optionally) additional annotation (see Details).
  }
  \item{sumprot}{
    Function (median or mean) for how the peptide values should be aggregated.
  }
}
\details{
This function is the main function called in \code{\link{customSummarizer}} and \code{\link{esetSummarizer}}.
The aggregate values for each variable in \code{exprs}, where the aggregation is determined by \code{sumprot}, are returned in a single data frame.

The argument \code{fdata} must contain a variable for protein IDs (called ``Protein''); 
    additional columns (e.g. sequence or peptide IDs) are allowed but are ignored.
}
\value{
  Data frame of the following variables:
  \describe{
    \item{\code{Protein}}{Factor; unique protein identifiers.}
    \item{\code{...    }}{   One variable for each of the original \code{exprs} data frame; aggregate values.}
    }
}
\seealso{
\code{\link{customSummarizer}}, \code{\link{esetSummarizer}}, \code{\link{pvalueSummarizer}}
}
\examples{
data(vanHoof)
sumdataVH <- summarizer(dataVH[,-c(1:2)], dataVH[,1:2])
head(sumdataVH)
}
\keyword{ methods }
