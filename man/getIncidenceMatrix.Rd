\name{getIncidenceMatrix}
\docType{methods}

\alias{getIncidenceMatrix}
\alias{getIncidenceMatrix,GeneSetCollection,character-method}
\alias{getIncidenceMatrix,GeneSetCollection,missing-method}
\alias{getIncidenceMatrix,DeconGeneSetCollection,ANY-method}
\alias{getIncidenceMatrix,list,character-method}
\alias{getIncidenceMatrix,list,missing-method}

\title{Making an incidence matrix for a gene set}

\description{An incidence matrix summarizes membership of gene identifiers
in gene sets.}

\usage{
  
  getIncidenceMatrix(x, features)
  
}

\arguments{
  \item{x}{An object of class \code{\link{GeneSetCollection}}, 
    \code{\link{DeconGeneSetCollection}}, or \code{list}, for which 
    an incidence matrix will be produced.}
  \item{features}{A optional character vector providing set of features for 
    which the incidence matrix will be produced. By default, the matrix will 
    include one column for every gene in the collection.}
}

\value{
  A matrix with as many rows as gene sets in the GeneSetCollection, and columns 
  for genes. Each value will be 1 where a gene is in set, otherwise 0. By
  default, there will be as many columns as there are unique genes across all
  sets in the collection. The user may specify a set of features for which an
  incidence matrix will be made. In that case, there will be as many columns as
  there are provided features.
}

\author{Jason A. Hackney}

\seealso{
    \code{\link{incidence}}
}

\examples{
    gsc <- DeconGeneSetCollection()
    
    ### Make an incidence matrix for the gene sets
    incMat <- getIncidenceMatrix(gsc)
    
    ### Make an incidence matrix for the gene sets,
    ### explicitly calling out which features to use.
    ### This is identical to the previous incidence matrix.
    allGenes <- unlist(unique(geneIds(gsc)))
    incMat <- getIncidenceMatrix(gsc, features = allGenes)
}