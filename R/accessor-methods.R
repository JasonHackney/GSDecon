## Accessor function for getting per-gene set summaries of expression values
## Only returns values for gene sets with p-value < pvalueCutoff

eigengenes <- function(object, pvalueCutoff = NULL)
{
    if(is.null(pvalueCutoff))
        pvalueCutoff = pvalueCutoff(object)
        
    pvs <- object@pvalues
    ret <- object@eigengenes[, pvs < pvalueCutoff, drop=FALSE]
    ret
}

pvalueCutoff <- function(object)
{
    return(object@pvalueCutoff)
}

"pvalueCutoff<-" <- function(object, value)
{
    object@pvalueCutoff <- value
    object@nComp <- sum(pvalues(object) < value)
    object
}

nComp <- function(object)
{
    return(object@nComp)
}

"nComp<-"  <- function(object, value)
{
    object@nComp <- value
    object
}

## What were the p-values for gene sets in the analysis
## that were below the pvalueCutoff

pvalues <- function(object)
{
    return(object@pvalues[object@pvalues < pvalueCutoff(object)])
}
