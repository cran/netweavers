\name{filterNetwork}
\alias{filterNetwork}
\title{
Filter a Network Based on Identifier
}
\description{
  The function \code{filterNetwork} filters a node-node format network to include only those 
  unique interactions that contain specified proteins and their first interactors. 
}
\usage{
filterNetwork(network, identifiers)
}
\arguments{
  \item{network}{
    Data frame or matrix; network in node-node format.
  }
  \item{identifiers}{
    Vector of identifiers.
  }
}
\details{
The first two columns of \code{network} should be protein or gene names (symbols, identifiers) in node-node network format; all other columns are ignored. 
The \code{identifiers} should be a vector with entries matching at least a subset of those in \code{network}.
The \code{network} is filtered to include only those unique rows with values in one of the first two columns matching
those of \code{identifiers}.
}
\value{
  Data frame.
}
\seealso{
\code{\link{runNetweavers}}
}
\examples{
\dontrun{
set.seed(1234)
net <- data.frame(intA = letters[1:10], intB = letters[11:20])
ids <- sample(letters[1:20],5)
filterNetwork(net, ids)
}
}
\keyword{ methods }
