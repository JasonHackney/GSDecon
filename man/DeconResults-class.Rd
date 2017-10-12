\name{DeconResults-class}
\docType{class}
\alias{DeconResults-class}
\alias{DeconResults}
\alias{nComp}
\alias{pvalueCutoff}
\alias{pvalueCutoff<-}
\alias{pvalues}

\title{Class ``DeconResults''}
\description{
  This class represents the results of running gene set-based deconvolution on a
  data set.
}

\section{Slots}{
  \describe{
    \item{\code{pvalueCutoff}:}{
      A single numeric value giving the significance cutoff at which gene sets
      are considered informative.
    }
    \item{\code{pvalues}:}{
      A numeric vector with an entry for each gene set that has a p-value less
      than \code{pvalueCutoff}
    }
    \item{\code{eigengenes}:}{
      A numeric matrix with one column for each significant gene set and one row
      for each sample in the expression data provided
    }
    \item{\code{nComp}:}{
      A numeric vector of length 1, giving the number of significant components
    }
  }
}

\author{J.A. Hackney}

\keyword{classes}
