\name{weightCalc}
\alias{weightCalc}
\title{
  Calculate Weights for NetWeAvers Algorithm
}
\description{
  Calculates a weight for each protein in an interaction network based on the
  number of neighbors the protein has in the network. 
}
\usage{
weightCalc(filt_net)
}
\arguments{
  \item{filt_net}{
    Data frame or a matrix containing at least two character, factor or numeric variables/columns representing a network in node-node format.
  }
}
\details{
  The argument \code{filt_net} should have the first two columns in node-node network format; all other columns are ignored. 
  The weight for a given protein is the inverse of the number of neighbors/interactors the protein has in the network.  
}
\value{
  Data frame on the following 2 variables:
  \describe{
    \item{\code{Protein}}{Factor or vector, depending on input, of protein IDs.}
    \item{\code{weight }}{Numeric vector of weights.}
  }
}
\seealso{
\code{\link{runNetweavers}}
}
\examples{
data(vanHoof)
wc <- weightCalc(filtNetVH)
head(wc)
}
\keyword{ methods }
