\name{deconComponents}
\docType{methods}

\alias{deconComponents}
\alias{deconComponents,ANY,GeneSet-method}
\alias{deconComponents,ANY,character-method}
\alias{deconComponents,ANY,NULL-method}
\alias{deconComponents,ANY,missing-method}
\alias{deconComponents,ANY,numeric-method}
\alias{deconComponents,EList,logical-method}
\alias{deconComponents,ExpressionSet,logical-method}

\title{Deriving an eigengene for a gene set given a gene expression matrix}

\description{
  This is the main function used by the \code{\link{decon}} method to estimate 
  relative cellularity. In general, this function is not meant to be called 
  directly. 
}

\usage{
    deconComponents(object, geneSet, nComp = 1, 
        method = c("SVD", "NMF", "ICA", "MDS"))
}

\arguments{
  \item{object}{An object of class \code{matrix}, \code{ExpressionSet}, or
   \code{EList}, specifying the sample expression values. For \code{matrix}
   objects, each column should represent a sample, each row a feature.}
  \item{geneSet}{An object of class \code{\link{GeneSet}},
    \code{\link{character}}, \code{\link{numeric}}. character values are
    assumed to be valid rownames of the expression matrix, as are geneIds of a
    GeneSet. Numeric values should be as long as the number of rows in the
    expression matrix, with a 1 where the gene is present in the gene set,
    otherwise 0. If it's NULL or missing, then all rows of the expression 
    matrix will be used.}
  \item{nComp}{Numeric of length 1, specifying how many components to
    calculate for the matrix. See details for discussion of how components are
    calculated.}
  \item{method}{Method used to estimate eigengenes. At present only "SVD" is 
    supported. Other methods should only be used for testing purposes.}
}

\details{ 
    This function calculates the gene set-wise relative expression values for
    the expression data found in \code{object}. This is typically used as a
    helper function for the \code{\link{decon}} function, which calls this on
    a per-gene set basis.

    Expression summaries are calculated as follows. First, genes are
    standardized to have zero mean and unit standard deviation, retaining the
    values used to center and scale the data. The singular value decomposition
    (svd) of the matrix is then calculated. For each component that is
    requested (\code{1:nComp}), all other eigenvalues in the \code{d} matrix
    of the svd is set to 0, and the svd matrices are multiplied to yield a
    normalized expression matrix. The values in the normalized expression
    matrix are then decentered and descaled using the values from the original
    expression matrix. The relative expression values for the component is the
    per-sample mean of this matrix.
}
\value{
    A matrix with as many rows as there are samples in \code{object}, and as
    many columns as specified by \code{nComp}. NOTE: Dimensions are not
    dropped, even when only a single component is returned.
}

\author{J.A. Hackney}

\seealso{
  \code{\link{svd}}, \code{\link{decon}}
}

\examples{
    \dontrun{
        library(GEOquery)
        library(hgu133plus2.db)
        
        ### Get tissue/cell specific gene sets
        deconGSC <- DeconGeneSetCollection()
        
        ### Load a dataset of mixed immune cell lines
        GSE11058 <- getGEO("GSE11058")[[1]]
        exprs(GSE11058) <- log2(exprs(GSE11058))
        annotation(GSE11058) <- "hgu133plus2"
        
        ### Map identifiers
        deconU133GSC <- mapIdentifiers(deconGSC, AnnotationIdentifier(),
            revmap(hgu133plus2ENTREZID))
        
        ### Just get the estimates for a single gene set, with none of the 
        ### bells and whistles of the decon function.
        tcellEigengene <- deconComponents(GSE11058, deconU133GSC[["T"]])
    }
}

