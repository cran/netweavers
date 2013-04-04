\name{findDenseClusters}
\alias{findDenseClusters}
\title{
  Find Dense Clusters in a Network
}
\description{
  This function first converts a data frame containing the node-node network format into an instance of the \code{igraph} class
  and then removes self loops and directionality. 
  Then it uses random walks to find dense, highly connected subgraphs, called clusters, which are returned in a data frame. 
}
\usage{
findDenseClusters(filt_net, identifiers, min_clus_size = 4, steps = 10)
}
\arguments{
  \item{filt_net}{
    Data frame or matrix of two columns; node-node network format.
  }
  \item{identifiers}{
    Numeric or character vector representing proteins used to filter \code{filt_net}.
  }
  \item{min_clus_size}{
    Integer indicating the required minimal size (number of vertices) of a cluster.  
  } 
  \item{steps}{
    Integer specifying the length of random walks to perform. 
  }
}
\details{
This function is meant for use in the NetWeAvers algorithm:
the first item in the output list can be used as the \code{clust_2_protein} argument in \code{\link{runNetweavers}}. 
The other items are for visualization using \pkg{igraph}.

The \code{identifiers} should be a vector with entries matching at least a subset of those in \code{filt_net}.

See \href{http://igraph.sourceforge.net/doc/R/walktrap.community.html}{\code{walktrap.community}} in the
\pkg{igraph} package for more details on \code{steps}, 
the parameter required for community-finding using random walks.
}
\value{
    A list of: 
    \item{clust_2_protein}{Data frame of the following variables (cluster information):
    \describe{
      \item{\code{clusID}}{  Factor; a cluster ID.}
      \item{\code{clusSIZE}}{Numeric vector of integers containing the size (number of proteins) in each cluster.}
      \item{\code{X1 }}{      Character vector of proteins; the first protein for each cluster.}
      \item{\code{...}}{     Additional cluster symbols, as many as \code{clusSIZE} dictates per cluster.}
    }}
    \item{subGraph}{A graph object with the edges deleted.}
    \item{denseClus}{The dense cluster membership in the form of a \href{http://igraph.sourceforge.net/doc/R/communities.html}{communities} object.}
}
\references{
\href{http://igraph.sourceforge.net/}{http://igraph.sourceforge.net/} 
}
\seealso{
\code{\link{runNetweavers}}, \code{\link{filterNetwork}}
}
\examples{
\dontrun{
data(vanHoof)
dc <- findDenseClusters(filtNetVH, pvalsSumVH$Protein)
str(dc, list.len=5)
}
}
\keyword{ methods }
\keyword{ models }
