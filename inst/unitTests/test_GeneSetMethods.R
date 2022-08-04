test_1DeconGeneSetCollection <- function()
{
    deconGSC <- DeconGeneSetCollection()
    checkTrue(class(deconGSC)=="DeconGeneSetCollection")
    checkTrue(length(deconGSC) == 30)
}

test_2IncidenceMatrices <- function()
{
    deconGSC <- DeconGeneSetCollection()
    
    incMat1 <- incidence(deconGSC)
    checkTrue(is(incMat1, "matrix"))
    checkTrue(nrow(incMat1) == 30)
    checkTrue(ncol(incMat1) == 1506)
    
    
    fts <- c("340267", "5178", "151306", "25900", "94", "3003", 
        "5962", "2958", "4628", "94274")
        
    incMat2 <- incidence(deconGSC, 
        features = fts)
    checkTrue(is(incMat2, "matrix"))
    checkTrue(nrow(incMat2) == 30)
    checkTrue(ncol(incMat2) == 10)
    checkTrue(all.equal(colnames(incMat2), fts))

    GSC <- GeneSetCollection(unlist(deconGSC))
    incMat3 <- getIncidenceMatrix(GSC, 
        features = fts)
    checkTrue(is(incMat3, "matrix"))
    checkTrue(nrow(incMat3) == 30)
    checkTrue(ncol(incMat3) == 10)
    checkTrue(all.equal(colnames(incMat3), fts))
}
