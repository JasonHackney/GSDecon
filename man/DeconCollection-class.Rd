\name{DeconCollection-class}
\docType{class}
\alias{DeconCollection-class}
\alias{DeconCollection}
\alias{GeneSetCollection,ExpressionSet,missing,DeconCollection-method}
\alias{GeneSetCollection,character,AnnotationIdentifier,DeconCollection-method}
\alias{GeneSetCollection,AnnotationIdentifier,DeconCollection-method}
\alias{GeneSetCollection,missing,AnnotationIdentifier,DeconCollection-method}
\alias{GeneSetCollection,missing,EntrezIdentifier,DeconCollection-method}
\alias{GeneSetCollection,ExpressionSet,AnnotationIdentifier,DeconCollection-method}

\title{Class ``DeconCollection''}
\description{
  This class is used for instantiation and manipulation of gene sets used for deconvolution-based analyses. These include immune cell types and normal expression from 35 different human tissues.
}

\section{Objects from the Class}{
Objects should be created using \code{DeconCollection}.
}

\section{Slots}{
  \describe{
    \item{\code{type}:}{
      Holds an object of class \code{ScalarCharacter} containing the
      string describing the collection type. In this case, \code{Decon}.
    }
    \item{\code{ids}:}{
      Nothing is currently done with any IDs stored here. Saved for 
      future use.
    }
  }
}

\value{
  An object of class \code{DeconCollection}, initialized with the appropriate type.
}

\author{J.A. Hackney}

\examples{
  dc <- DeconCollection()
  deconGSC <- GeneSetCollection(idType=EntrezIdentifier(), setType=dc)
}
\keyword{classes}
