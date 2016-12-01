# Methods for the SlalomModel class

#' Create a new SlalomModel object.
#'
#' Slalom fits relatively complicated hierarchical Bayesian factor analysis
#' models with data and results stored in a \code{"SlalomModel"} object. This
#' function builds a new \code{"SlalomModel"} object from minimal inputs.
#'
#' @param object \code{"SCESet"} object N x G expression data matrix (cells x genes)
#' @param genesets a \code{"GeneSetCollection"} object containing annotated
#' gene sets
#' @param n_param number of hidden factors to fit in the model (2-5 recommended)
#' @param prune_genes logical, should genes that are not annotated to any gene
#' sets be filtered out?
#' @param min_genes scalar, minimum number of genes required in order to retain
#' a gene set for analysis

#' @return a new SlalomModel object
#'
#' @details
#' This function builds and returns the object, checking for validity, which
#' includes checking that the input data is of consistent dimensions.
#'
#' @import GSEABase
#' @import scater
#' @importFrom Biobase exprs
#' @importFrom methods new
#' @importFrom methods validObject
#' @export
#' @examples
#' gmtfile <- system.file("extdata", "reactome_subset.gmt", package = "slalom")
#' genesets <- GSEABase::getGmt(gmtfile)
#' data("mesc")
#' model <- newSlalomModel(mesc, genesets, n_hidden = 5, min_genes = 10)
newSlalomModel <- function(object, genesets, n_hidden = 5, prune_genes = TRUE,
                           min_genes = 15) {
    Y <- t(Biobase::exprs(object))
    ## convert GeneSetCollection into I matrix of indicators
    I <- lapply(genesets, function(x) {colnames(Y) %in% GSEABase::geneIds(x)})
    names(I) <- names(genesets)
    I <- as.matrix(as.data.frame(I))
    rownames(I) <- colnames(Y)
    ## filter genesets on min_genes
    n_sets_pass <- colSums(I) >= min_genes
    if (sum(n_sets_pass) > 0L)
        I <- I[, n_sets_pass]
    else
        stop("No gene sets found containing more than min_genes genes.")
    ## prune genes if desired
    if (prune_genes) {
        keep_gene <- (rowSums(I) > 0L)
        I <- I[keep_gene, ]
    }
    retained_genes <- rownames(I)
    ## add hidden factors
    if (n_hidden > 0L) {
        hidden_factors <- matrix(TRUE, nrow = nrow(I), ncol = n_hidden)
        colnames(hidden_factors) <- paste0("hidden", sprintf("%02d", 1:n_hidden))
        I <- cbind(I, hidden_factors)
    }

    ## create new SlalomModel object for output
    out <- new("SlalomModel",
        K = ncol(I),
        N = nrow(Y),
        G = nrow(I),
        Y = Y[, retained_genes],
        I = I,
        termNames = colnames(I),
        geneNames = retained_genes,
        cellNames = scater::cellNames(object),
        iUnannotatedDense = c(rep(FALSE, ncol(I) - n_hidden),
                              rep(TRUE, n_hidden))
    )

    ## Check validity of object
    validObject(out)
    out
}


# genevars <- matrixStats::colVars(mesc)
# keepgene <- genevars > sort(genevars)[3000]
# sceset <- newSCESet(exprsData = t(mesc[, keepgene]))
# mesc <- sceset
# save(mesc, file = "../../data/mesc.RData")
# object <- mesc
#
# gl <- lapply(g, function(x) {geneIds(x) <- capwords(geneIds(x), strict = TRUE); x})
# gl <- GeneSetCollection(gl)
# gsmall <- gl[grep("CELL", names(gl))]
# toGmt(gsmall, con = "reactome_subset.gmt")
# genesets <- gsmall


################################################################################
### Define validity check for SlalomModel class object

setValidity("SlalomModel", function(object) {
    msg <- NULL
    valid <- TRUE

    if (valid) TRUE else msg
})