# tests for SlalomModel methods

test_that("newSlalomModel produces a valid object", {
    gmtfile <- system.file("extdata", "reactome_subset.gmt", package = "slalom")
    genesets <- GSEABase::getGmt(gmtfile)
    data("mesc")
    model <- newSlalomModel(mesc, genesets, n_hidden = 5, min_genes = 10)
})
