
setClass("DeconCollection",
    contains = "CollectionIdType",
    representation = representation(
        type = "ScalarCharacter"
    ),
    prototype = prototype(
        type = new("ScalarCharacter", "Decon")
    )
)

DeconCollection <- function()
{
    new("DeconCollection")
}

setClass("DeconGeneSetCollection",
    contains = "GeneSetCollection"
)

setClass("DeconResults",
    representation = representation(
        pvalueCutoff = "numeric",
        pvalues = "numeric",
        eigengenes = "matrix",
        nComp = "numeric"
    )
)

setClass("DeconGeneSet",
    contains = "GeneSet",
    representation = representation(
        geneValues = "numeric"
    ),
    validity = function(object)
    {
        res <- FALSE
        ## geneValues must be the same length as geneIds
        if(length(geneIds(object)) == length(geneValues(object))) res <- TRUE
        res
    }
)
