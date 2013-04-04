\name{fillup}
\alias{fillup}
\title{
  Combine Vectors of Varying Length into a Data Frame 
}
\description{
  Combines elements of a list, each of which are vectors with possibly 
  varying lengths, into a data frame. 
}
\usage{
fillup(mylist)
}
\arguments{
  \item{mylist}{
    A list.
  }
}
\details{
  This function takes a list of vectors with possibly varying lengths and 
  returns a single data frame where each row is one of the vectors from \code{mylist}
  and missing values are filled in with \code{NA}. The row order is defined by the 
  order in \code{mylist}.
}
\value{
  A data frame containing aggregated values from \code{mylist}.   
}
\seealso{
\code{\link{findDenseClusters}}
}
\examples{
x <- list(2:6,1:10,5:7)
fillup(x)
}
\keyword{ methods }
