### Take in expression data and generate summarized values for a gene set
### The resulting data have one value per sample per component requested
### The default method is used mostly for the signature 
### object = "numeric", geneSet = "logical", with most other cases handled
### below

setGeneric("deconComponents",
    signature = signature("object", "geneSet"),
    def = function(object, geneSet, nComp = 1, 
        method = c("SVD", "NMF", "ICA", "MDS"))
    {
        method = match.arg(method, c("SVD", "NMF", "ICA", "MDS"))
        
        if(nComp > 0){
            if(nrow(object[geneSet,,drop=FALSE]) < 2){
                svs <- matrix(NA, nrow = ncol(object), ncol = 1)
                warning("A gene set with too few genes was used. Gene sets must have >=2 genes. Returning NA.")
            } else {
    
                if(method == "SVD"){
                    object <- t(scale(t(object[geneSet,])))
                    scl <- attributes(object)$"scaled:scale"
                    cnt <- attributes(object)$"scaled:center"
                
                    ## Use SVD to compute the fatorization of the expression
                    ## matrix and pull in the particular eigenvectors that are
                    ## interesting from the data.
                    svd.res <- svd(object)
                    svs <- matrix(0, nrow = ncol(object), ncol = nComp)
                    for(i in 1:nComp){
                        svs[, i] <- eigencomponent(svd.res, i, 
                            center = cnt, scale = scl)
                    }
            
                } else if(method == "ICA"){
                    object <- t(scale(t(object[geneSet,])))
        
                    ## Use ICA to compute a factorization of the 
                    ## expression matrix
                    ica.res <- fastICA(object, n.comp = length(nComp), 
                        method = "C")
                    # signs <- vapply(colSums(ica.res$S > 0)/nrow(ica.res$S), 
                    #     function(x) ifelse(x > 0.5, 1, -1), numeric(1))
                    # svs <- t(t(ica.res$K) * signs)
                    svs <- ica.res$K
        
                } else if(method == "NMF"){
                    object <- object[geneSet,]

                    ## Use non-negative matrix factorization on the matrix to
                    ## factorize the matrix into the appropriate number of
                    ## factors
                    nmf.res <- nmf(object, rank = nComp)
                    svs <- t(.coef(nmf.res))
        
                }else if(method == "MDS"){
                    object <- t(scale(t(object[geneSet,])))
                    # dmat <- dist(object)
                    dmat <- as.dist(1-cor(object))
                    ## Use MDS to compute the fatorization of the expression
                    ## matrix and pull in the particular eigenvectors that are
                    ## interesting from the data.
                    svs <- cmdscale(dmat, k = nComp)
                }
            } 
        } else {
            svs <- NULL
        }
        svs
    }
)

### GeneSet method: if no genes are in common between the GeneSet and
### expression data, throw an error.

setMethod("deconComponents",
    signature = signature(
        object = "ANY",
        geneSet = "GeneSet"
    ),
    definition = function(object, geneSet, nComp, method = "SVD")
    {
        geneSet <- intersect(rownames(object), geneIds(geneSet))
        if(length(geneSet) == 0){
            stop("No genes in the gene set were found as rownames of the expression data.")
        }
        res <- deconComponents(object, geneSet, nComp, method = method)
        res
    }
)

### character method: if no genes are in common between the character vector and
### expression data, throw an error.

setMethod("deconComponents",
    signature = signature(
        object = "ANY",
        geneSet = "character"
    ),
    definition = function(object, geneSet, nComp, method = "SVD")
    {
        geneSet <- intersect(rownames(object), geneSet)
        if(length(geneSet) == 0){
            stop("No genes in the gene set were found as rownames of the expression data.")
        }
        geneSet <- rownames(object) %in% geneSet
        res <- deconComponents(object, geneSet, nComp, method = method)
        res
    }
)

### If the gene set is NULL, then just use the whole matrix
### hopefully the user has already subsetted the expression data...

setMethod("deconComponents",
    signature = signature(
        object = "ANY",
        geneSet = "NULL"
    ),
    definition = function(object, geneSet, nComp, method = "SVD")
    {
        geneSet <- rep(TRUE, nrow(object))
        res <- deconComponents(object, geneSet, nComp, method = method)
        res
    }
)

### If the gene set is missing, then just use the whole matrix
### hopefully the user has already subsetted the expression data...

setMethod("deconComponents",
    signature = signature(
        object = "ANY",
        geneSet = "missing"
    ),
    definition = function(object, geneSet, nComp, method = "SVD")
    {
        geneSet <- rep(TRUE, nrow(object))
        res <- deconComponents(object, geneSet, nComp, method = method)
        res
    }
)

### If the gene set is numeric, assume that what we got was an incidence
### matrix, with all 0's and 1's. This would mean geneSet is the same length
### as the number of rows in the expression data. If not, maybe it was a list
### of indices into the expression data. But check to make sure we don't have
### an overrun here.

setMethod("deconComponents",
    signature = signature(
        object = "ANY",
        geneSet = "numeric"
    ),
    definition = function(object, geneSet, nComp, method = "SVD")
    {
    
        if(all(geneSet %in% c(0,1))){
            ### First guess: it's just an incidence matrix: all values are 
            ### 0 or 1 and the length is the same as the number of rows in 
            ### the expression data.
            if(length(geneSet) == nrow(object)){
                geneSet <- geneSet == 1
            } else {
                stop("Incidence matrix dimensions don't match expression data.")
            }
        } else {
            ### Second guess, they are indices of rows in the matrix. 
            if(! all(geneSet %in% 1:nrow(object))){
                ### But maybe not.
                stop("Got indices too large for number of features in expression data.")
            } else {
                ### If so, make a logical vector out of it
                incMat <- rep(FALSE, nrow(object))
                incMat[geneSet] <- TRUE
                geneSet <- incMat
            }
        }
    
        res <- deconComponents(object, geneSet, nComp, method = method)
        res
    }
)


setMethod("deconComponents",
    signature = signature(
        object = "ExpressionSet",
        geneSet = "logical"
    ),
    definition = function(object, geneSet, nComp, method = "SVD")
    {
        object <- exprs(object)
        res <- deconComponents(object, geneSet, nComp, method = method)
        res
    }
)


setMethod("deconComponents",
    signature = signature(
        object = "EList",
        geneSet = "logical"
    ),
    definition = function(object, geneSet, nComp, method = "SVD")
    {
        object <- object$E
        res <- deconComponents(object, geneSet, nComp, method = method)
        res
    }
)


rescale <- function(x)
{
    rng <- range(x, na.rm=TRUE)
    x <- (x - rng[1] + 0.01)/(rng[2] - rng[1] + 0.01)
    x
}

eigencomponent <- function(svd, n = 1, center = 0, scale = 1)
{
    newD <- svd$d
    newD[-n] <- 0
    sv <- svd$u %*% diag(newD) %*% t(svd$v)
    sv <- sweep(sweep(sv, 1, FUN = "*", scale), 1, FUN = "+", center)
    sv <- colMeans(sv)
    sv
}
