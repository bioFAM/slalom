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
    rdsfile <- system.file("extdata", "sim_N_20_v3.rds", package = "slalom")
    sim <- readRDS(rdsfile)
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(logcounts = sim[["init"]][["Y"]])
    )
    rownames(sce) <- unique(unlist(GSEABase::geneIds(genesets[1:20])))[1:500]
    m <- newSlalomModel(sce, genesets[1:20], n_hidden = 1, min_genes = 1)
    ## initialise this model
    m <- init(m, pi_prior = sim[["init"]][["Pi"]], n_hidden = 1)
    m$X_E1 <- sim[["init"]][["X"]]
    m$W_E1 <- sim[["init"]][["W"]]
    m$alpha_E1 <- sim[["init"]][["Alpha"]]
    m$epsilon_E1 <- sim[["init"]][["Epsilon"]]
    m$W_gamma0 <- sim[["init"]][["W_gamma0"]]

    expect_that(m, is_a("Rcpp_SlalomModel"))

    ## test precise results from one update
    m <- updateSlalom(m)
    expect_equal(m$alpha_E1[,1], sim[["first_iter"]][["Alpha"]])
    expect_equal(m$X_E1, sim[["first_iter"]][["X"]])
    expect_equal(m$W_E1, sim[["first_iter"]][["W"]])
    expect_equal(m$epsilon_E1[,1], sim[["first_iter"]][["Epsilon"]])
    expect_equal(m$W_gamma0, sim[["first_iter"]][["W_gamma0"]])

    ## train for 500 iterations and compare to Python results
    m <- train(m, minIterations = 400, nIterations = 500, shuffle = FALSE)
    expect_true(cor(c(m$alpha_E1[,1]),
                    c(sim[["final_iter"]][["Alpha"]])) > 0.9999)
    expect_true(cor(c(m$X_E1), c(sim[["final_iter"]][["X"]])) > 0.999999)
    expect_true(cor(c(m$W_E1), c(sim[["final_iter"]][["W"]])) > 0.999999)
    expect_true(cor(c(m$epsilon_E1[,1]),
                    c(sim[["final_iter"]][["Epsilon"]])) > 0.999999)
    expect_true(cor(c(m$W_gamma0), c(sim[["final_iter"]][["W_gamma0"]])) >
                    0.999999)

    # par(mfrow = c(3, 7), mar = c(0.1, 0.1, 0.1, 0.1))
    # for (i in 1:21) {
    #     plot(m$W_E1[,i], sim[["final_iter"]][["W"]][,i])
    #     text(-1, 0, labels = i)
    #     abline(0, 1, col = "firebrick")
    # }
    # for (i in 1:21) {
    #     plot(m$W_gamma0[,i], sim[["final_iter"]][["W_gamma0"]][,i])
    #     text(-1, 0, labels = i)
    #     abline(0, 1, col = "firebrick")
    # }
    # for (i in 1:21) {
    #     plot(m$X_E1[,i], sim[["final_iter"]][["X"]][,i])
    #     text(-1, 0, labels = i)
    #     abline(0, 1, col = "firebrick")
    # }
    # par(mfcol = c(1, 2))
    # plot(m$alpha_E1[,1], sim[["final_iter"]][["Alpha"]])
    # abline(0, 1, col = "firebrick")
    # plot(m$epsilon_E1[,1], sim[["final_iter"]][["Epsilon"]])
    # abline(0, 1, col = "firebrick")

    ## access Python results if needed
    # library(rhdf5)
    # h5file <- file.path("./inst", "extdata", "outSim", "N_20_20it_v3.h5py")
    # h5ls(h5file)



    ###########################################################################
    gmtfile <- system.file("extdata", "c2.cp.reactome.v5.2.symbols.gmt",
                           package = "slalom")
    #gmtfile <- file.path("./inst", "extdata", "c2.cp.reactome.v5.2.symbols.gmt")
    genesets <- GSEABase::getGmt(gmtfile)
    genenames <- unique(unlist(GSEABase::geneIds(genesets)))
    ## mesc data
    data("mesc")
    rownames(mesc) <- toupper(rownames(mesc))
    table(rownames(mesc) %in% genenames)
    model <- newSlalomModel(mesc, genesets[1:100], n_hidden = 5, min_genes = 10)
    model <- init(model, seed = 222)
    expect_that(model, is_a("Rcpp_SlalomModel"))

    model <- train(model, nIterations = 10)

    expect_that(model, is_a("Rcpp_SlalomModel"))

})


