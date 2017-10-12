\name{eigengenes}
\docType{methods}

\alias{eigengenes}

\title{Retrieving eigengenes from DeconResults}

\description{
  This is the accessor function for retrieving the eigengenes calculated by
  the decon method. 
}

\usage{
    eigengenes(object, pvalueCutoff = NULL)
}

\arguments{
  \item{object}{An object of class \code{\link{DeconResults}}.}
  \item{pvalueCutoff}{numeric value giving the p-value cutoff, below 
    which eigegenes are returned. If \code{NULL}, then the \code{pvalueCutoff} 
    from the DeconResults object is used.}
}
\value{
    A matrix with as many rows as there are samples in \code{object}, and as
    many columns as there were gene sets provided to \code{\link{decon}}.
}

\author{J.A. Hackney}

\seealso{
  \code{\link{DeconResults}}, \code{\link{decon}}
}

\examples{
    \dontrun{
        library(GEOquery)
        library(hgu133plus2.db)
        
        ### Get tissue/cell specific gene sets
        deconGSC <- DeconGeneSetCollection()
        
        ### Dataset with mixed immune cell lines
        GSE11058 <- getGEO("GSE11058")[[1]]
        exprs(GSE11058) <- log2(exprs(GSE11058))
        annotation(GSE11058) <- "hgu133plus2"
        
        ### Map the identifiers to Affymetrix probe sets
        deconU133GSC <- mapIdentifiers(deconGSC, AnnotationIdentifier(),
            revmap(hgu133plus2ENTREZID))
        
        ### Run the deconvolution algorithm
        deconResults <- decon(GSE11058, ~1, deconU133GSC)
        
        ### Retrieve the results
        eigengenes <- eigengenes(deconResults, pvalueCutoff = 0.05)
    }
}