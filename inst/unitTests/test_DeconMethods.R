deconGSC <- DeconGeneSetCollection()
library(GEOquery)
library(genefilter)
GSE11058 <- getGEO("GSE11058")[[1]]
exprs(GSE11058) <- log2(exprs(GSE11058))
annotation(GSE11058) <- "hgu133plus2"
GSE11058 <- featureFilter(GSE11058)

cellFractions <- cbind(
    Jurkat  = c(1,1,1,0,0,0,0,0,0,0,0,0,0.25,0.25,0.25,0.05,0.05,0.05,0.01,0.01,
                0.01,0.002,0.002,0.002),
    "IM-9"  = c(0,0,0,1,1,1,0,0,0,0,0,0,0.125,0.125,0.125,0.317,0.317,0.317,
                0.495,0.495,0.495,0.333,0.333,0.333),
    Raji    = c(0,0,0,0,0,0,1,1,1,0,0,0,0.25,0.25,0.25,0.475,0.475,0.475,0.165,
                0.165,0.165,0.333,0.333,0.333),
    "THP-1" = c(0,0,0,0,0,0,0,0,0,1,1,1,0.375,0.375,0.375,0.158,0.158,0.158,
                0.33,0.33,0.33,0.333,0.333,0.333),
    "B cell" = rep(NA, 24))

rownames(cellFractions) <- c("GSM279589", "GSM279590", "GSM279591",
     "GSM279592", "GSM279593", "GSM279594", "GSM279595", "GSM279596",
     "GSM279597", "GSM279598", "GSM279599", "GSM279600", "GSM279601",
     "GSM279602", "GSM279603", "GSM279604", "GSM279605", "GSM279606",
     "GSM279607", "GSM279608", "GSM279609", "GSM279610", "GSM279611",
     "GSM279612")
cellFractions[, "B cell"] <- rowSums(
    cellFractions[, 2:3])/rowSums(cellFractions[, 1:4])

GSE11058 <- GSE11058[, rownames(cellFractions)]
deconU133GSC <- mapIdentifiers(deconGSC, AnnotationIdentifier(),
    revmap(hgu133plus2ENTREZID))

test_3DeconMethods <- function()
{
    GSE11058_model <- model.matrix(~1, GSE11058)
    
    ### Test full decon method with all three required arguments
    deconRes1 <- decon(GSE11058, GSE11058_model, deconU133GSC, doPerm = FALSE)
    checkTrue(class(deconRes1) == "DeconResults")
    
    ### Test full decon method with formula
    deconRes2 <- decon(GSE11058, ~1, geneSets = deconU133GSC, doPerm = FALSE)
    checkTrue(class(deconRes2) == "DeconResults")
    
    ### Test full decon method with incidence matrix
    incMat <- incidence(deconU133GSC, features = rownames(GSE11058))
    deconRes3 <- decon(GSE11058, ~1, incMat, doPerm = FALSE)
    checkTrue(class(deconRes3) == "DeconResults")
    
    ### Test full decon method with one gene set
    deconRes4 <- decon(GSE11058, GSE11058_model, deconU133GSC[1], doPerm = FALSE)
    checkTrue(class(deconRes4) == "DeconResults")
    
    eigs4 <- eigengenes(deconRes4)
    checkTrue(length(dim(eigs4)) == 2)
    checkTrue(ncol(eigs4) == 1)
    
    ### Test full decon method with missing design
    deconRes5 <- decon(GSE11058, geneSets = deconU133GSC, doPerm = FALSE)
    checkTrue(class(deconRes5) == "DeconResults")
    
    ### Test full decon method with missing design and incidence matrix
    deconRes5 <- decon(GSE11058, geneSets = incMat, doPerm = FALSE)
    checkTrue(class(deconRes5) == "DeconResults")
    
    ### Check dimensions of eigengenes
    eigs <- eigengenes(deconRes1)
    checkTrue(nrow(eigs) == ncol(GSE11058))
    checkTrue(ncol(eigs) == length(deconU133GSC))
    
}

test_4DeconComponentsMethods <- function()
{
    fts <- geneIds(deconU133GSC[["T"]])
    
    ### ExpressionSet, geneSet method
    vals1 <- deconComponents(GSE11058, geneIds(deconU133GSC[["T"]]))
    checkTrue(length(vals1) == 24)
    
    ### matrix, geneSet method
    vals2 <- deconComponents(exprs(GSE11058), geneIds(deconU133GSC[["T"]]))
    checkTrue(length(vals2) == 24)
    
    ### matrix, character method
    vals3 <- deconComponents(exprs(GSE11058), fts)
    checkTrue(length(vals3) == 24)
    
    ### matrix, logical method
    vals4 <- deconComponents(exprs(GSE11058), rownames(GSE11058) %in% fts)
    checkTrue(length(vals4) == 24)
    
    ### matrix, numeric method
    incMat <- incidence(deconU133GSC, rownames(GSE11058))
    vals5 <- deconComponents(exprs(GSE11058), incMat["T", ])
    checkTrue(length(vals5) == 24)
    
    allVals <- cbind(vals1, vals2, vals3, vals4, vals5)
    nuniq <- apply(allVals, 1, function(x) length(unique(x) == 1))
    checkTrue(all(nuniq == 1))

}