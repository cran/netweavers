\name{vanHoof}
\alias{vanHoof}
\alias{dataVH}
\alias{phenoDataVH}
\alias{ExpSetVH}
\alias{summarizedExpSetVH}
\alias{DEtestVH}
\alias{proteinInfoVH}
\alias{filtNetVH}
\alias{clust_2_proteinVH}
\alias{networkVH}
\alias{pvalsSumVH}
\docType{data}
\title{
Van Hoof Dataset
}
\description{
Example dataset from Van Hoof \emph{et al.} (2009), including raw data (data frame and 
\href{http://bioconductor.org/packages/release/bioc/vignettes/Biobase/inst/doc/ExpressionSetIntroduction.pdf}{\code{ExpressionSet}}), 
phenotype data, summarized data (\code{ExpressionSet}), protein \emph{p}-values and other statistics used in NetWeAvers, 
linear model statistics, a filtered protein-protein interaction network, and dense clusters of proteins.
}
\usage{
data(vanHoof)
}
\format{
  There are 10 objects in the \sQuote{vanHoof} dataset.
  \describe{
    \item{\code{dataVH}}{            Data frame of the mass spectrometry data, pre-processed using PVIEW (Khan \emph{et al.} 2009, Khan \emph{et al.} 2011).} 
    \item{\code{phenoDataVH}}{       Data frame of the experimental design for \code{dataVH}.}
    \item{\code{ExpSetVH}}{          An \code{ExpressionSet} containing \code{dataVH} and \code{phenoDataVH}.}
    \item{\code{summarizedExpSetVH}}{An \code{ExpressionSet}; the output of \code{\link{esetSummarizer}} using \code{ExpSetVH}.}
    \item{\code{DEtestVH}}{          List; the output of \code{\link{DEtest}} using \code{ExpSetVH}.}
    \item{\code{networkVH}}{         Data frame; Reactome version 43 protein-protein interaction network in node-node format 
    (\href{http://www.reactome.org/download/current/homo_sapiens.interactions.txt.gz}{http://www.reactome.org/download/current/homo_sapiens.interactions.txt.gz)}.}
	\item{\code{filtNetVH}}{         Data frame; the \code{\link{networkVH}} filtered based on the proteins in \code{dataVH}.}
    \item{\code{clust_2_proteinVH}}{ Data frame; the first element of the list output by \code{\link{runNetweavers}} on \code{filtNetVH} and \code{pvalsSumVH}.}  
    \item{\code{proteinInfoVH}}{     Data frame; the second element of the list output by \code{\link{runNetweavers}} on \code{filtNetVH} and \code{pvalsSumVH}.}
	\item{\code{pvalsSumVH}}{        Summarized \emph{p}-values from \code{DEtestVH} comparison of 30 to 60 minutes.}
  }
}
\source{
The authors of Van Hoof \emph{et al.} 2009.
}
\references{
Khan Z, Bloom JS, Garcia BA, Singh M, and Kruglyak L. (2009). Protein Quantification Across Hundreds of Experimental Conditions. \emph{Proc. Natl. Acad. Sci.} 06(37):15544-15548.

Khan Z, Amini S, Bloom JS, Ruse C, Caudy AA, Kruglyak L, Singh M, Perlman DH, and Tavazoie S. (2011). Accurate proteome-wide protein quantification from high-resolution 15N mass spectra. \emph{Genome Biology} 12:R122.

Van Hoof D, Munoz J, Braam SR, Pinkse MW, Linding R, Heck AJ, Mummery CL, and Krijgsveld J. (2009). Phosphorylation dynamics during early differentiation of human embryonic stem cells. \emph{Cell Stem Cell} 5(2):214-226.
}
\examples{
data(vanHoof)
}
\keyword{datasets}
