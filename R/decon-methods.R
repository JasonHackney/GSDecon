setGeneric("decon",
    def = function(object, model = NULL, geneSets, doPerm = TRUE, nPerm = 249,
        pvalueCutoff = 0.01, nComp = 1, trim = FALSE, seed = NULL, ...)
    {
        standardGeneric("decon")
    },
    signature = c("object", "model", "geneSets")
)

### If you get a DeconGeneSetCollection, make it into an incidence matrix
### At least some of the identifiers in the gene sets should be present in the 
### expression data, though not all need to be

setMethod("decon",
    signature = signature(object = "ANY", model = "ANY", 
        geneSets = "DeconGeneSetCollection"),
    definition = function(object, model = NULL, geneSets, doPerm = TRUE, nPerm = 249,
        pvalueCutoff = 0.01, nComp = 1, trim = FALSE, seed = NULL, ...)
    {
        ugenes <- unlist(geneIds(geneSets))
        
        if(!any(ugenes %in% rownames(object))){
            stop("Gene sets do not have any features in common with expression data. 
Check that the annotations are the same between your gene sets and expression data.")
        }
        
        gs.mat <- incidence(geneSets, rownames(object))

        res <- decon(object, model, geneSets = gs.mat, 
            doPerm = doPerm, nPerm = nPerm, pvalueCutoff = pvalueCutoff, 
            nComp = nComp, trim = trim, seed = seed, ...)

        res
    }
)

### If you get a list for geneSets, assume it's a list of identifiers, and 
### make an incidence matrix of it.

setMethod("decon",
    signature = signature(object = "ANY", model = "ANY", 
        geneSets = "list"),
    definition = function(object, model = NULL, geneSets, doPerm = TRUE, nPerm = 249,
        pvalueCutoff = 0.01, nComp = 1, trim = FALSE, seed = NULL, ...)
    {
        ugenes <- unlist(geneSets)
        
        if(!any(ugenes %in% rownames(object))){
            stop("Gene sets do not have any features in common with expression data. 
Check that the annotations are the same between your gene sets and expression data.")
        }
        
        gs.mat <- getIncidenceMatrix(geneSets, rownames(object))
            res <- decon(object, model, geneSets = gs.mat, 
            doPerm = doPerm, nPerm = nPerm, pvalueCutoff = pvalueCutoff, 
            nComp = nComp, trim = trim, seed = seed, ...)
        res
    }
)


## For eSet-based objects, giving a formula for the model, we should be able
## to make that into a matrix using model.matrix

setMethod("decon",
    signature = signature(object = "eSet", model = "formula",
        geneSets = "matrix"),
    definition = function(object, model = NULL, geneSets, doPerm = TRUE, nPerm = 249,
        pvalueCutoff = 0.01, nComp = 1, trim = FALSE, seed = NULL, ...)
    {
        mm <- model.matrix(model, pData(object))
        res <- decon(object, mm, geneSets, 
            doPerm = doPerm, nPerm = nPerm, pvalueCutoff = pvalueCutoff, 
            nComp = nComp, trim = trim, seed = seed, ...)
        res
    }
)

## If we get all matrices, assume that you can just use limma to run a linear
## model on the specified design matrix. Most of the work is done by
## deconValuesAndComponents

setMethod("decon",
    signature = signature(object = "matrix", model = "matrix", 
        geneSets = "matrix"),
    definition = function(object, model = NULL, geneSets, doPerm = TRUE, nPerm = 249,
        pvalueCutoff = 0.01, nComp = 1, trim = FALSE, seed = NULL)
    {
        if(sum(geneSets) == 0){
            stop("Gene sets do not have any features in common with expression data. 
Check that the annotations are the same between your gene sets and expression data.")
        }
        ret <- deconValuesAndComponents(object, model, geneSets, 
            doPerm = doPerm, nPerm = nPerm, pvalueCutoff = pvalueCutoff, 
            nComp = nComp, trim = trim, seed = seed)
        return(ret)
    }
)


## ExpressionSet methods -- The easiest ones of all. Just take use the 
## exprs matrix to compute the different components

setMethod("decon",
    signature = signature(object = "ExpressionSet", model = "missing", 
        geneSets = "matrix"),
    definition = function(object, model = NULL, geneSets, doPerm = TRUE, nPerm = 249,
        pvalueCutoff = 0.01, nComp = 1, trim = FALSE, seed = NULL, ...)
    {    
        model <- model.matrix(~1, object)
        ret <- decon(object, model, geneSets, 
            doPerm = doPerm, nPerm = nPerm, pvalueCutoff = pvalueCutoff, 
            nComp = nComp, trim = trim, seed = seed, ...)
        return(ret)
    }
)

setMethod("decon",
    signature = signature(object = "ExpressionSet", model = "NULL", 
        geneSets = "matrix"),
    definition = function(object, model = NULL, geneSets, doPerm = TRUE, nPerm = 249,
        pvalueCutoff = 0.01, nComp = 1, trim = FALSE, seed = NULL, ...)
    {    
        model <- model.matrix(~1, object)
        ret <- decon(object, model, geneSets, 
            doPerm = doPerm, nPerm = nPerm, pvalueCutoff = pvalueCutoff, 
            nComp = nComp, trim = trim, seed = seed, ...)
        return(ret)
    }
)

setMethod("decon",
    signature = signature(object = "ExpressionSet", model = "matrix", 
        geneSets = "matrix"),
    definition = function(object, model = NULL, geneSets, doPerm = TRUE, nPerm = 249,
        pvalueCutoff = 0.01, nComp = 1, trim = FALSE, seed = NULL, ...)
    {    
        if(sum(geneSets) == 0){
            stop("Gene sets do not have any features in common with expression data. 
Check that the annotations are the same between your gene sets and expression data.")
        }
        object <- exprs(object)
        ret <- deconValuesAndComponents(object, model, geneSets, 
            doPerm = doPerm, nPerm = nPerm, pvalueCutoff = pvalueCutoff, 
            nComp = nComp, trim = trim, seed = seed)
        return(ret)
    }
)


## Count based methods for CountDataSet, DESeqDataSet and DGEList. 
## Use voom to get variance stabilized data, then fit a model using limma to 
## get the residuals. The variance stabilized data will also be used for 
## calculating the surrogate variables

setMethod("decon",
    signature(object = "CountDataSet", model = "missing", 
        geneSets = "matrix"),
    definition = function(object, model = NULL, geneSets, doPerm = TRUE, nPerm = 249,
        pvalueCutoff = 0.01, nComp = 1, trim = FALSE, seed = NULL, ...)
    {
        ## Need to find the design of the CountDataSet
        ## If missing, assume just using conditions from the object
        mm <- model.matrix(~condition, pData(object))
        res <- decon(object, mm, geneSets, 
            doPerm = doPerm, nPerm = nPerm, pvalueCutoff = pvalueCutoff, 
            nComp = nComp, trim = trim, seed = seed, ...)
        res
    }
)

setMethod("decon",
    signature(object = "CountDataSet", model = "matrix", 
        geneSets = "matrix"),
    definition = function(object, model = NULL, geneSets, doPerm = TRUE, nPerm = 249,
        pvalueCutoff = 0.01, nComp = 1, trim = FALSE, seed = NULL, ...)
    {
        ## Just use the model.matrix as-is
        dge <- DGEList(counts(object))
        res <- decon(dge, model, geneSets, 
            doPerm = doPerm, nPerm = nPerm, pvalueCutoff = pvalueCutoff, 
            nComp = nComp, trim = trim, seed = seed, ...)
        res
    }
)

setMethod("decon",
    signature(object = "DESeqDataSet", model = "missing", 
        geneSets = "matrix"),
    definition = function(object, model = NULL, geneSets, doPerm = TRUE, nPerm = 249,
        pvalueCutoff = 0.01, nComp = 1, trim = FALSE, seed = NULL, ...)
    {
        mm <- model.matrix(design(object), as.data.frame(colData(object)))
        dge <- DGEList(counts(object))
        res <- decon(dge, mm, geneSets, 
            doPerm = doPerm, nPerm = nPerm, pvalueCutoff = pvalueCutoff, 
            nComp = nComp, trim = trim, seed = seed, ...)
        res
    }
)

setMethod("decon",
    signature(object = "DESeqDataSet", model = "matrix", 
        geneSets = "matrix"),
    definition = function(object, model = NULL, geneSets, doPerm = TRUE, nPerm = 249,
        pvalueCutoff = 0.01, nComp = 1, trim = FALSE, seed = NULL, ...)
    {
        dge <- DGEList(counts(object))
        res <- decon(dge, model, geneSets, 
            doPerm = doPerm, nPerm = nPerm, pvalueCutoff = pvalueCutoff, 
            nComp = nComp, trim = trim, seed = seed, ...)
        res
    }
)

setMethod("decon",
    signature(object = "DESeqDataSet", model = "formula", 
        geneSets = "matrix"),
    definition = function(object, model = NULL, geneSets, doPerm = TRUE, nPerm = 249,
        pvalueCutoff = 0.01, nComp = 1, trim = FALSE, seed = NULL, ...)
    {
        mm <- model.matrix(model, as.data.frame(colData(object)))
        dge <- DGEList(counts(object))
        res <- decon(dge, mm, geneSets, 
            doPerm = doPerm, nPerm = nPerm, pvalueCutoff = pvalueCutoff, 
            nComp = nComp, trim = trim, seed = seed, ...)
        res
    }
)


setMethod("decon",
    signature(object = "DGEList", model = "matrix", geneSets = "matrix"),
    definition = function(object, model = NULL, geneSets, doPerm = TRUE, nPerm = 249,
        pvalueCutoff = 0.01, nComp = 1, trim = FALSE, seed = NULL)
    {
        if(sum(geneSets) == 0){
            stop("Gene sets do not have any features in common with expression data. 
Check that the annotations are the same between your gene sets and expression data.")
        }
        ## Just use the model.matrix as-is
        object <- calcNormFactors(object)
        vsd <- voom(object, design = model)
        res <- deconValuesAndComponents(vsd, model, geneSets, doPerm = doPerm, 
            nPerm = nPerm, pvalueCutoff = pvalueCutoff, nComp = nComp, 
            trim = trim, seed = seed)
    }
)

setMethod("decon",
    signature(object = "EList", model = "missing", geneSets = "matrix"),
    definition = function(object, model = NULL, geneSets, doPerm = TRUE, nPerm = 249,
        pvalueCutoff = 0.01, nComp = 1, trim = FALSE, seed = NULL)
    {
        if(sum(geneSets) == 0){
            stop("Gene sets do not have any features in common with expression data. 
Check that the annotations are the same between your gene sets and expression data.")
        }
        rsds <- apply(object$E, 1, sd, na.rm=TRUE)
        if(any(rsds == 0)){
            stop("Genes with zero variance in the input dataset. Please remove these and try again.")
        }
        
        ## If no model matrix is present, use the design from the EList
        model <- object$design
        res <- deconValuesAndComponents(object, model, geneSets, 
            doPerm = doPerm, nPerm = nPerm, pvalueCutoff = pvalueCutoff, 
            nComp = nComp, trim = trim, seed = seed)
    }
)

setMethod("decon",
    signature(object = "EList", model = "matrix", geneSets = "matrix"),
    definition = function(object, model = NULL, geneSets, doPerm = TRUE, nPerm = 249,
        pvalueCutoff = 0.01, nComp = 1, trim = FALSE, seed = NULL)
    {
        if(sum(geneSets) == 0){
            stop("Gene sets do not have any features in common with expression data. 
Check that the annotations are the same between your gene sets and expression data.")
        }
        
        rsds <- apply(object$E, 1, sd, na.rm=TRUE)
        if(any(rsds == 0)){
            stop("Genes with zero variance in the input dataset. Please remove these and try again.")
        }

        ## Just use the model.matrix as-is
        res <- deconValuesAndComponents(object, model, geneSets, 
            doPerm = doPerm, nPerm = nPerm, pvalueCutoff = pvalueCutoff, 
            nComp = nComp, trim = trim, seed = seed)
    }
)

### Summarize the expression data for gene sets:
### First do permutation testing on the residuals matrix from fitting the null 
### model. Then calculate the scores for each gene set using deconComponents

deconValuesAndComponents <- function(object, model, geneSets, doPerm = TRUE, 
    nPerm = 249, pvalueCutoff = 0.01, nComp = 1, trim = FALSE, seed = NULL,
    method = c("SVD", "NMF", "ICA", "MDS"))
{
    if (!is.null(seed)) {
        set.seed(seed)
    }
    
    if(nrow(object) != ncol(geneSets)){
        stop("Gene set incidence matrix and expression data have different numbers of features.")
    }
    
    if(ncol(object) != nrow(model)){
        stop("Expression data and design matrix have different numbers of samples")
    }
    
    method = match.arg(method, c("SVD", "NMF", "ICA", "MDS"))
    
    if(class(object) == "EList"){
        rsds <- rowSds(object$E)
    } else {
        rsds <- rowSds(object)
    }
    if(any(rsds == 0)){
        object <- object[rsds > 0,]
        warning("Genes with 0 variance removed.")
    }
    
    if(trim){
        geneSets <- t(apply(geneSets, 1, trimGeneSet, dataSet = object))
    }
               
    resid <- residuals(lmFit(object, model), object)
    nComp <- rep(nComp, nrow(geneSets))
    nComp <- pmin(nComp, rowSums(geneSets == 1))
    nm <- make.names(rep(rownames(geneSets), nComp), unique=TRUE)
    
    ## Center and scale the residuals, so I can save some
    ## time on larger datasets
    res.pvs <- lapply(1:nrow(geneSets), function(x) 
        permuteEigenvalues(geneSets[x,], resid=resid, 
        doPerm=doPerm, nPerm=nPerm, nComp = nComp[x], 
        gsName = rownames(geneSets)[x]))
    res.pvs <- do.call("c", res.pvs)
    names(res.pvs) <- nm
    
    gs.svs <- lapply(1:nrow(geneSets), function(x) deconComponents(object, 
        geneSets[x, ], method = method, nComp = nComp[x]))
    gs.svs <- do.call("cbind", gs.svs)
    colnames(gs.svs) <- nm
    rownames(gs.svs) <- colnames(object)
    
    nComp <- sum(res.pvs < pvalueCutoff)
    
    retval <- new("DeconResults", pvalues = res.pvs, 
        pvalueCutoff = pvalueCutoff, eigengenes = gs.svs, nComp = nComp)
    retval

}

# =============================================================================
# Calculate the first n eigenvalues of a subset of a matrix, and then do 
# resampling to take other same-sized subsets of the matrix and calculate their 
# first n eigenvalues. The p-value is the number of observed p-values greater 
# than that observed for the real subset.
# -----------------------------------------------------------------------------

permuteEigenvalues <- function(resid, gs, doPerm = TRUE, nPerm = 249, 
    seed = NULL, nComp = 1, pvalueCutoff = 0.01, gsName = NULL)
{
    if(sum(gs == 1) < 2){
        retval <- rep(NA, nComp)
    } else { 
        # If someone says they want to do 0 permutations, let them
        if(nPerm == 0){
            nPerm <- 100
            doPerm <- FALSE
        }
        
        gs.svd <- svd(resid[gs == 1, ])
        ## Calculate the eigenvalues for the original matrix
        origVal <- gs.svd$d[1:(nComp+2)]
        pvs <- rep(1/(nPerm+1), nComp+2)
        
        if(doPerm){
            if(!missing(seed) & !is.null(seed))
                set.seed(seed)
            #Permute membership in the gene set
            perm.vals <- do.call("cbind",
                lapply(c(1:nPerm), function(x) {
                    perm.svd <- list(d=NA)
                    try(
                        perm.svd <- svd(resid[sample(gs, 
                            replace=FALSE) == 1, ]),
                        silent=TRUE
                    )
                    ## Get the permuted eigenvalues
                    retval <- perm.svd$d[1:(nComp+2)]
                    retval
                })
            )
            if(any(is.na(perm.vals))){
                warning("Something has gone wrong in the permutations. \
                    Not all permutations were completed.")}
            pvs <- ((rowSums(origVal <= perm.vals, na.rm=TRUE)+1) / 
                (rowSums(!is.na(perm.vals))+1))

            if(sum(pvs < pvalueCutoff) > nComp){
                warning(paste("More than", nComp, "eigenvalues were significant in", gsName))
            }
        }
        retval <- pvs[1:nComp]
    }
    retval
}
