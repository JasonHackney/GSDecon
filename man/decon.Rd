\name{decon}
\docType{methods}

\alias{decon}
\alias{decon,ANY,ANY,DeconGeneSetCollection-method}
\alias{decon,ANY,ANY,list-method}
\alias{decon,CountDataSet,formula,matrix-method}
\alias{decon,CountDataSet,missing,matrix-method}
\alias{decon,CountDataSet,matrix,matrix-method}
\alias{decon,DESeqDataSet,formula,matrix-method}
\alias{decon,DESeqDataSet,missing,matrix-method}
\alias{decon,DESeqDataSet,matrix,matrix-method}
\alias{decon,DGEList,matrix,matrix-method}
\alias{decon,EList,matrix,matrix-method}
\alias{decon,EList,missing,matrix-method}
\alias{decon,ExpressionSet,matrix,matrix-method}
\alias{decon,ExpressionSet,missing,matrix-method}
\alias{decon,ExpressionSet,NULL,matrix-method}
\alias{decon,eSet,formula,matrix-method}
\alias{decon,matrix,matrix,matrix-method}

\title{Identification of tissue-specific surrogate variables using pre-specified
  gene sets}

\description{
  Estimate surrogate variables in a data matrix from cell type- or tissue-
  specific gene sets. The surrogate variables are determined by the
  \code{\link{deconComponents}} method.
}

\usage{
  decon(object, model = NULL, geneSets, doPerm = TRUE, nPerm = 249, 
      pvalueCutoff = 0.01, nComp = 1, trim = FALSE, seed = NULL, ...)  
}

\arguments{
  \item{object}{An object of class \code{matrix}, \code{ExpressionSet},
    \code{CountDataSet}, \code{DESeqDataSet}, or \code{DGEList}, specifying the
    sample expression values. For \code{matrix} objects, each column should
    represent a sample, each row a feature, and it is assumed that using the
    \code{lmFit} function from limma would be appropriate for the data in the 
    matrix. Count data are first transformed to log2-scale data using the 
    \code{\link{voom}} function from the \code{limma} package}
  \item{model}{A formula or an n x k design matrix specifying the model of
    interest see \code{\link{model.matrix}}. Default value is NULL. For objects     that are \code{CountDataSet} or \code{DESeqDataSet} object, in which case
    the design of the object is used to create a model matrix using the
    \code{\link{model.matrix}} function. If the object is an
    \code{ExpressionSet}, then a model matrix is created using an
    intercept-only model.}
  \item{geneSets}{An object of class \code{\link{DeconGeneSetCollection}}, a
    \code{list} of character vectors, or an incidence matrix with dimensions 
    g x m, where the columns represent genes, and the rows represent gene sets,
    with a 1 where a gene is in a gene set, otherwise 0. For \code{list}, the
    values should correspond to the row names of the expression data. For
    \code{\link{DeconGeneSetCollection}}, the \code{\link{geneIds}} should
    correspond to row names of the expression data.}
  \item{doPerm}{A boolean value of length 1, specifying if permutation testing
    for significance of gene sets should be performed. See below for
    details on the permutation testing and its interpretation}
  \item{nPerm}{A numeric value of length 1, specifying how many permutations 
    should be performed}
  \item{pvalueCutoff}{A numeric value of length 1, generally between 1/nPerm and
    1 specifying at what significance level should gene sets be considered
    informative in the dataset.}
  \item{nComp}{A numeric value of length 1, specifying how many components to
    test for in each gene set. For well-formed gene sets, this should be
    set to 1, indicating that the gene set should have one major set of
    correlated genes. If more components that nComp are found to be
    significant, a warning is thrown, as the gene set is likely not well
    specified.}
  \item{trim}{Logical. Should gene sets be trimmed before summarization. If
    TRUE, then only genes with an average pairwise correlation coefficient >
    0.1 are included in the decon algorithm. Otherwise, the whole gene set is
    used.}
  \item{seed}{A seed to set for random number generator used in the permutation.
    Setting this will allow for reproducible p-values to be generated for the 
    gene sets.}
  \item{\dots}{Currently not used, but may be used in the future.}
}

\details{ 
    \code{decon} attempts to identify gene sets that are significantly
    informative in the residuals matrix for a given linear model. This is done
    through a gene-wise permutation strategy. For each permutation, the first
    \code{nComp} eigenvalues are compared to the eigenvalues of the original
    gene set expression matrix. An empirical p-value is calculated by finding
    how many random eigenvalues are greater than the observed eigenvalues.

    For gene sets that are considered significant (permutation p-value less
    than the supplied alpha), an eigengene for the gene set is calculated. For
    each significant eigenvalue, an eigengene is calculated by the method
    described for \code{\link{deconComponents}}. The first eigengene typically
    represents the relative amount of that cell or tissue type in the mixed
    sample.

    Eigengenes beyond the first can be somewhat difficult to interpret, and by
    default are not looked for. However, looking for significance of the
    second (or third) eigenvalue can be informative about the relative
    consistency within a gene set. Ideally, the majority of the variance of
    the gene set would be explained by the first eigenvector. If there is a
    large amount of variance explained by the second eigenvector, this
    suggests that your gene set is identifying two separate expression
    patterns in the data set of interest.    
}

\value{
    An object of class \code{\link{DeconResults}} with the following slots:
    \item{pvalueCutoff}{A single numeric value giving the significance cutoff at
      which gene sets are considered informative}
    \item{pvalues}{A numeric vector with an entry for each gene set that has a
      p-value less than \code{pvalueCutoff}}
    \item{eigengenes}{A numeric matrix with one column for each significant gene 
      set and one row for each sample in the expression data provided}
    \item{nComp}{A numeric vector of length 1, giving the number of significant 
      components}
}

\author{J.A. Hackney}

\seealso{
  \code{\link{svd}}
  \code{\link{deconComponents}}
}

\examples{
    \dontrun{
        library(GEOquery)
        library(hgu133plus2.db)
        
        deconGSC <- DeconGeneSetCollection()
        GSE11058 <- getGEO("GSE11058")[[1]]
        exprs(GSE11058) <- log2(exprs(GSE11058))
        annotation(GSE11058) <- "hgu133plus2"
        
        deconU133GSC <- mapIdentifiers(deconGSC, AnnotationIdentifier(),
            revmap(hgu133plus2ENTREZID))
        deconResults <- decon(GSE11058, ~1, deconU133GSC)
    }
}
