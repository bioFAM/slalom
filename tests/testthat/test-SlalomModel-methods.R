# tests for SlalomModel methods

library(slalom)

test_that("newSlalomModel produces a valid object", {
    gmtfile <- system.file("extdata", "reactome_subset.gmt", package = "slalom")
    genesets <- GSEABase::getGmt(gmtfile)
    data("mesc")
    model <- newSlalomModel(mesc, genesets, n_hidden = 5, min_genes = 10)
    expect_that(model, is_a("Rcpp_SlalomModel"))
})


test_that("SlalomModel initialises", {
    gmtfile <- system.file("extdata", "reactome_subset.gmt", package = "slalom")
    genesets <- GSEABase::getGmt(gmtfile)
    data("mesc")
    model <- newSlalomModel(mesc, genesets, n_hidden = 5, min_genes = 10)
    model <- init(model)
    expect_that(model, is_a("Rcpp_SlalomModel"))

})


test_that("SlalomModel trains", {
    gmtfile <- system.file("extdata", "reactome_subset.gmt", package = "slalom")
    genesets <- GSEABase::getGmt(gmtfile)

    ###########################################################################
    ## simulated data
    rdsfile <- system.file("extdata", "sim_N_20_v2.rds", package = "slalom")
    sim <- readRDS(rdsfile)
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(logcounts = sim[["init"]][["Y"]])
    )
    rownames(sce) <- unique(unlist(GSEABase::geneIds(genesets[1:20])))[1:500]
    m <- newSlalomModel(sce, genesets[1:20], n_hidden = 1, min_genes = 1)
    m <- init(m, pi_prior = sim[["init"]][["Pi"]], n_hidden = 1)
    m$X_E1 <- sim[["init"]][["X"]]
    m$W_E1 <- sim[["init"]][["W"]]
    m$alpha_E1 <- sim[["init"]][["Alpha"]]
    m$epsilon_E1 <- sim[["init"]][["Epsilon"]]
    m$W_gamma0 <- sim[["init"]][["W_gamma0"]]

    ## test precise results from one update
    m <- updateSlalom(m)
    expect_equal(m$alpha_E1[,1], sim[["first_iter"]][["Alpha"]])
    expect_equal(m$X_E1, sim[["first_iter"]][["X"]])
    expect_equal(m$W_E1, sim[["first_iter"]][["W"]])
    expect_equal(m$epsilon_E1[,1], sim[["first_iter"]][["Epsilon"]])

    ### train this model
    m <- newSlalomModel(sce, genesets[1:20], n_hidden = 1, min_genes = 1)
    m <- init(m, pi_prior = sim[["init"]][["Pi"]], n_hidden = 1)
    m$X_E1 <- sim[["init"]][["X"]]
    m$W_E1 <- sim[["init"]][["W"]]
    m$alpha_E1 <- sim[["init"]][["Alpha"]]
    m$epsilon_E1 <- sim[["init"]][["Epsilon"]]
    m$W_gamma0 <- sim[["init"]][["W_gamma0"]]

    m <- train(m, nIterations = 1500)

    expect_that(m, is_a("Rcpp_SlalomModel"))

    ###########################################################################
    ## mesc data
    data("mesc")
    model <- newSlalomModel(mesc, genesets, n_hidden = 5, min_genes = 10)
    model <- init(model, seed = 222)
    expect_that(model, is_a("Rcpp_SlalomModel"))

    model <- train(model, nIterations = 5)

    expect_that(model, is_a("Rcpp_SlalomModel"))

})


