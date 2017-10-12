
# =============================================================================
# For convenience, just return the description of all the GeneSets in a
# GeneSetCollection
# -----------------------------------------------------------------------------

setMethod("description",
         signature=signature(object='DeconGeneSetCollection'),
         definition=function(object, ...)
    {
        res <- lapply(object, description)
        names(res) <- lapply(object, setName)
        return(res)
    }
)

setMethod("description<-",
         signature=signature(object='DeconGeneSetCollection',
            value='character'),
         definition=function(object, value)
    {
        desc <- function (x, y) 
        {
            slot(x, "shortDescription") <- mkScalar(y)
            `slot<-`(x, "setIdentifier", check = FALSE, 
                mkScalar(.getGUID()))
            validObject(x)
            x
        }
        objs <- GeneSetCollection(mapply(desc, object, value))
        objs
        
    }
)

setMethod("collectionType<-",
         signature=signature(object='DeconGeneSetCollection', value='CollectionType'),
         definition=function(object, value)
    {
        colType <- function (x, y) 
        {
            slot(x, "collectionType") <- y
            `slot<-`(x, "setIdentifier", check = FALSE,
                mkScalar(.getGUID()))
            validObject(x)
            x
        }
        value <- rep(list(value), length.out=length(object))
        objs <- GeneSetCollection(mapply(colType, object, value))
        objs

    }
)

.getGUID <- function(gs){
    guid <- digest(list(setName(gs), description(gs),
        gs@version, geneIdType(gs), collectionType(gs),
        organism(gs), creationDate(gs)))
    guid
}

# =============================================================================
# A new incidence function. Deprecate this in favor of the next one!
# -----------------------------------------------------------------------------

setGeneric('getIncidenceMatrix', def=function(x, features) 
    standardGeneric('getIncidenceMatrix'))

setMethod("getIncidenceMatrix",
    signature = signature(x="GeneSetCollection", features="character"),
    definition = function(x, features)
    {
        mat <- matrix(data=0, ncol=length(features),nrow=length(x))
        geneids <- geneIds(x)
        for(i in seq_along(x)){
            mat[i,match(geneids[[i]], features)] <- 1
        }
        colnames(mat) <- features
        rownames(mat) <- vapply(x, setName, character(1))
        mat
    }
)

setMethod("getIncidenceMatrix",
    signature = signature(x="GeneSetCollection", features="missing"),
    definition = function(x, features)
    {
        mat <- incidence(x)
        mat
    }
)

setMethod("getIncidenceMatrix",
    signature = signature(x = "DeconGeneSetCollection", features = "ANY"),
    definition = function(x, features)
    {
        mat <- incidence(x, features = features)
        mat
    }
)

setMethod("getIncidenceMatrix",
    signature = signature(x = "list", features = "character"),
    definition = function(x, features)
    {
        myMat <- matrix(0, ncol = length(features), nrow = length(x))
        colnames(myMat) <- features
        rownames(myMat) <- names(x)
        for(gsi in seq_along(x)){
            myMat[gsi, x[[gsi]] ] <- 1
        }
        myMat
    }
)

setMethod("getIncidenceMatrix",
    signature = signature(x = "list", features = "missing"),
    definition = function(x, features)
    {
        features <- unique(unlist(x))
        myMat <- getIncidenceMatrix(x, features = features)
        myMat
    }
)

# =============================================================================
# The incidence method for an DeconGeneSetCollection is like the original
# incidence method, but could return a larger number of columns, which
# might not be present as geneIds in the Collection
# -----------------------------------------------------------------------------

setMethod("incidence", 
    signature = signature(x = "DeconGeneSetCollection"),
    definition = function(x, features = unique(unlist(geneIds(x))))
    {
        mat <- matrix(data=0, ncol=length(features), nrow=length(x))
        geneids <- geneIds(x)
        for(i in seq_along(x)){
            mat[i,match(geneids[[i]], features)] <- 1
        }
        colnames(mat) <- features
        rownames(mat) <- vapply(x, setName, character(1))
        mat
    }
)

# =============================================================================
# For some reason, you can only map individual GeneSets using a bimap or 
# environment. I want to do more with mapping
# -----------------------------------------------------------------------------

.mapIds <- function(gs, idType, bimapList, features=NULL, ...)
{
    foundIds <- unique(unlist(bimapList[intersect(geneIds(gs), 
        names(bimapList))]))
    if(length(foundIds) > 0){
        foundIds <- foundIds[!is.na(foundIds)]
        if( ! is.null(features)) foundIds <- intersect(foundIds, features)
    }
    if( is.null(foundIds) || length(foundIds) == 0) foundIds <- character()
    # Make a new gene set here. Copy over the attributes, and give it the right 
    # idType
    new.gs <- gs
    geneIds(new.gs) <- foundIds
    new.gs@geneIdType <- idType
    new.gs
}

## Using the base mapIdentifiers from GSEABase is way too slow. I'll just use
## my own function, tyvm

setMethod('mapIdentifiers',
    signature = signature(what="DeconGeneSetCollection", 
        to="GeneIdentifierType", from="AnnDbBimap"),
    definition = function(what, to, from, ..., verbose = FALSE)
    {
        fromList <- mget(unique(unlist(geneIds(what))), from, ifnotfound=NA)
        gsc <- GeneSetCollection(lapply(what, .mapIds, to, fromList, ...))
        gsc <- new("DeconGeneSetCollection", gsc)
        gsc
    }
)

setMethod('mapIdentifiers',
    signature = signature(what="DeconGeneSetCollection", 
        to="GeneIdentifierType", from="environment"),
    definition = function(what, to, from, ..., verbose = FALSE)
    {
        gsc <- GeneSetCollection(lapply(what, mapIdentifiers, to, from, ..., 
            verbose = verbose))
        gsc <- new("DeconGeneSetCollection", gsc)
        gsc
    }
)

# =============================================================================
# How to instantiate an Decon GeneSetCollection
# -----------------------------------------------------------------------------

.getBaseGSC <- function(){
    gsc <- readRDS(system.file('extdata', 'decon_gsc.rds', package='GSDecon'))
}

.getMap <- function(idType)
{
    map <- revmap(getAnnMap('ENTREZID', annotation(idType)))
    map
}

setMethod('GeneSetCollection',
    signature = signature(
        object = 'missing',
        idType = 'EntrezIdentifier',
        setType = 'DeconCollection'
    ),
    function(object, ..., idType, setType)
    {
        decon.gsc <- .getBaseGSC()
        decon.gsc
    }
)

setMethod('GeneSetCollection',
    signature = signature(
        object = 'missing',
        idType = 'AnnotationIdentifier',
        setType = 'DeconCollection'
    ),
    function(object, ..., idType, setType)
    {
        decon.gsc <- .getBaseGSC()
        map <- .getMap(idType)
        gsc.mapped <- mapIdentifiers(decon.gsc, idType, map)
        gsc.mapped
    }
)

setMethod('GeneSetCollection',
    signature = signature(
        object = 'ExpressionSet',
        idType = 'missing',
        setType = 'DeconCollection'
    ),
    function(object, ..., idType, setType)
    {
        decon.gsc <- .getBaseGSC()
        idType <- AnnotationIdentifier(annotation(object))
        map <- .getMap(idType)
        gsc.mapped <- mapIdentifiers(decon.gsc, idType, map,
            features=featureNames(object))
        gsc.mapped
    }
)

setMethod('GeneSetCollection',
    signature = signature(
        object = 'ExpressionSet',
        idType = 'AnnotationIdentifier',
        setType = 'DeconCollection'
    ),
    function(object, ..., idType, setType)
    {
        decon.gsc <- .getBaseGSC()
        map <- .getMap(idType)
        gsc.mapped <- mapIdentifiers(decon.gsc, idType, map,
            features=featureNames(object))
        gsc.mapped
    }
)

setMethod('GeneSetCollection',
    signature = signature(
        object = 'character',
        idType = 'AnnotationIdentifier',
        setType = 'DeconCollection'
    ),
    function(object, ..., idType, setType)
    {
        decon.gsc <- .getBaseGSC()
        map <- .getMap(idType)
        gsc.mapped <- mapIdentifiers(decon.gsc, idType, map,
            features=object)
        gsc.mapped
    }
)


