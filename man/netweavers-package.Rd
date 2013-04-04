\name{netweavers-package}
\alias{netweavers-package}
\alias{netweavers}
\docType{package}
\title{
NetWeAvers: Weighted Averages for Networks
}
\description{
This package finds dense clusters of protein-protein interactions (PPIs) based on 
\emph{p}-values from proteomics data.
Peptide values may be summarized or rolled-up to protein values, 
a differential expression test is provided to obtain \emph{p}-values for each protein, 
and these \emph{p}-values, weighted by 
number of interactors in a network, may be averaged to get 
scores for dense clusters in a network.
A permutation test is provided to assess the significance level of a cluster.
}
\details{
The function \code{\link{runNetweavers}} runs the NetWeAvers algorithm as a wrapper for all of the main functions in the package.
Each step of the NetWeAvers algorithm is offered as a separate function so they can be run independently. This is useful, for example, 
when a user wants to vary parameter values, would like to decrease computation time by avoiding repeating certain steps or 
if only the output from one of the steps (e.g. the differential test) is desired.   

Two objects (a PPI network and protein \emph{p}-values) are required for the algorithm, 
along with the specification of a few parameters, including whether or not to run a
permutation test.

The PPI network should be in node-node format (e.g. proteinA proteinB). A dense cluster-finding algorithm
segregates the network into a few smaller subnetworks that have high connectivity. This part of NetWeAvers can be run separately using 
\code{\link{findDenseClusters}}, the output of which can be input into \code{\link{runNetweavers}} instead of the PPI network itself.

The protein \emph{p}-values for NetWeAvers can be generated from any testing procedure, 
but this package includes a differential expression test (\code{\link{DEtest}}) for quantified peptide- or protein-level data. At
the protein level, the output from the function may be directly used in \code{\link{runNetweavers}}.
Two options are available for peptide-level data: 
\enumerate{
  \item{First test for differential expression on the peptide-level using \code{\link{DEtest}} 
        and then get a summarized \emph{p}-value for each protein using Fisher's method (\code{\link{pvalueSummarizer}}).} 
  \item{Summarize or roll-up the peptide-level data to protein-level data using \code{\link{customSummarizer}} or \code{\link{esetSummarizer}} 
       (depending on whether your data are in matrix or ExpressionSet form) and then test at the protein-level using \code{\link{DEtest}}.}
}
The output from either option can be directly input in the NetWeAvers algorithm.
}
\author{
Elizabeth A. McClellan

Maintainer: Elizabeth A. McClellan <e.mcclellan@erasmusmc.nl>
}
\keyword{ package }
