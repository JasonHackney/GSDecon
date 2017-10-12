### Methods for making DeconGeneSetCollections.
### We want this class essentially just to get some 
### specialized use out of the incidence method and 
### some nicer mapping of identifiers

setGeneric('DeconGeneSetCollection', 
    def=function(object, idType = EntrezIdentifier(), 
        setType = DeconCollection())
    {
        standardGeneric('DeconGeneSetCollection')
    },
    signature = signature("object")
)

### If we got a list, see if we can make that a GeneSetCollection
### If so, recast to a DeconGeneSetCollection

setMethod("DeconGeneSetCollection",
    signature = signature(object = "list"),
    definition = function(object, idType, setType)
    {
        gsc <- GeneSetCollection(object)
        decon.gsc <- new("DeconGeneSetCollection", gsc)
        decon.gsc
    }
)

### If we got a GeneSetCollection, just need to recast to 
### DeconGeneSetCollection.

setMethod("DeconGeneSetCollection",
    signature = signature(object = "GeneSetCollection"),
    definition = function(object, idType, setType)
    {
        decon.gsc <- new("DeconGeneSetCollection", object)
        decon.gsc
    }
)

### If we didn't get something, but got ENtrezIdentifier and 
### DeconCollection as a type, then we can just pull the default
### gene sets included with the package and make those into a 
### DeconGSC

setMethod("DeconGeneSetCollection",
    signature = signature(object = "missing"),
    definition = function(object, idType = EntrezIdentifier(), 
        setType = DeconCollection())
    {
        gsc <- GeneSetCollection(idType = idType, setType = setType)
        decon.gsc <- new("DeconGeneSetCollection", gsc)
        decon.gsc
    }
)
    
### If we get anything else, try to make a GeneSetCollection
### out of it, and see if that works.

setMethod("DeconGeneSetCollection",
    signature = signature(object = "ANY"),
    definition = function(object, idType, setType)
    {
        gsc <- GeneSetCollection(object, idType = idType, setType = setType)
        decon.gsc <- new("DeconGeneSetCollection", gsc)
        decon.gsc
    }
)

