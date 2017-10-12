\name{DeconGeneSetCollection}
\docType{class}

\alias{DeconGeneSetCollection}
\alias{description}
\alias{description<-}
\alias{collectionType<-}
\alias{show}

\alias{DeconGeneSetCollection}
\alias{DeconGeneSetCollection-class}
\alias{DeconGeneSetCollection,ANY-method}
\alias{DeconGeneSetCollection,GeneSetCollection-method}
\alias{DeconGeneSetCollection,list-method}
\alias{DeconGeneSetCollection,missing-method}

\alias{description,DeconGeneSetCollection-method}
\alias{description<-,DeconGeneSetCollection,character-method}

\alias{collectionType<-,DeconGeneSetCollection,CollectionType-method}

\alias{show,DeconGeneSetCollection-method}

\title{Collection of cell- or tissue-specific gene sets}

\description{
  Classes used for collections of cell- or tissue-specific gene sets for use in deconvolution, along with specialized accessor functions. The primary use for this class is the use of a specialized incidence function, and some slightly improved mapping capabilities.
}

\usage{
    DeconGeneSetCollection(object, idType = EntrezIdentifier(), 
        setType = DeconCollection())
    
    description(object, ...)
    description(object) <- value
    
    collectionType(object) <- value

}

\arguments{
  \item{object}{A \code{\link{GeneSetCollection}}, or any other object that can be coerced to a \code{\link{GeneSetCollection}} using the constructor for that class.}
  \item{idType}{An argument of class `GeneIdentifierType', used to indicate how the `geneIds' will be represented.}
  \item{setType}{An argument of class `CollectionType', used to indicate what kind of collection is created. Defaults to \code{IRISCollection}, the default class for cell type specific gene sets.}
  \item{value}{ For \code{description}, a character vector as long as the number of \code{GeneSet}s in the collection, for \code{collectionType}, a  vector of \code{\link{CollectionType}}s as long as the number of \code{GeneSet}s in the collection}
  \item{\dots}{Further arguments to be passed on to other methods.}
}

\value{
  For \code{DeconGeneSetCollection}, an object of class \code{DeconGeneSetCollection} with geneIds from the appropriate geneIdType. The gene sets in the collection represent genes that are relatively specific to that particular tissue or cell type.
  For \code{description}, a character vector with the description from each of the \code{GeneSet}s in the collection.
}

\seealso{
  \code{\link{DeconCollection-class}}, \code{\link{GeneSetCollection-class}}
}

\examples{
    
    ## Get a set of tissue specific marker genes.
    dgsc <- DeconGeneSetCollection()
    
    ### Get the descriptions from the collection
    description(dgsc)
    
    ### Make a normal GeneSetCollection
    gs1 <- GeneSet(setName = "Set1", shortDescription = "Example 1")
    gs2 <- GeneSet(setName = "Set2", shortDescription = "Example 2")
    gsc2 <- GeneSetCollection(list(gs1, gs2))
    ### Coerce it to a DeconGeneSetCollection
    dgsc2 <- DeconGeneSetCollection(gsc2)
    dgsc2
    
    ### Map them to a popular microarray platform
    library(hgu133plus2.db)
    u133GSC <- mapIdentifiers(dgsc, AnnotationIdentifier(), hgu133plus2ENTREZID)
}