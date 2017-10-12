
setMethod("show",
    signature = signature(object = "DeconGeneSetCollection"),
    definition = function(object)
        {
            some <- function(x) paste(paste(Biobase::selectSome(x, 4),
                collapse = ", "), " (", length(x), " total)", sep = "")
            gids <- unique(unlist(geneIds(object)))
            itypes <- unique(vapply(lapply(object, geneIdType), class, character(1)))
            ctypes <- unique(vapply(lapply(object, collectionType), class, character(1)))
            cat("DeconGeneSetCollection\n", "  names: ", some(names(object)),
                "\n", "  unique identifiers: ", some(gids), "\n", 
                "  types in collection:\n", "    geneIdType: ", some(itypes), 
                "\n", "    collectionType: ", some(ctypes), "\n", sep = "")
        }
)


setMethod("show",
    signature = signature(object = "DeconResults"),
    definition = function(object)
    {
        pvs <- pvalues(object)
        selpvs <- names(pvs[pvs < pvalueCutoff(object)])
        
        cat("DeconResults\n", 
            "  pvalueCutoff: ", pvalueCutoff(object), "\n",
            "  significant components: ", 
                paste(Biobase::selectSome(selpvs), collapse = ", "), 
                " (", nComp(object), " total)\n",
            sep = ""
        )
    }
)
