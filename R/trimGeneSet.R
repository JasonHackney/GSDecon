### Given a gene set, which genes should we include?
### Just atke those that have some minimal level of pair-wise correlation
### across the entire dataset.

trimGeneSet <- function(geneSet, dataSet, cor.cutoff = 0.1)
{
    if (sum(geneSet == 1) == 0) { stop("no genes in input gene set") }
    if (cor.cutoff < -1 | cor.cutoff > 1) { 
        stop("wrong input for correlation") 
    }
    
    if(is(dataSet, "ExpressionSet")){
        dataSet <- exprs(dataSet)
    } else if( is(dataSet, "EList") ){
        dataSet <- dataSet$E
    }
        
    corTab <- cor(t(dataSet[geneSet == 1, ]), method="spearman")
    diag(corTab) <- NA
    meanCor <- rowMeans(corTab, na.rm=TRUE)
    
    geneSet[geneSet == 1][meanCor < cor.cutoff] <- 0
    geneSet
}
