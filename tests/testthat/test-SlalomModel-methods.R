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

    library(rhdf5)
    smallh5 <- system.file("extdata", "outSim", "G_500N_20",
                         "data.h5py", package = "slalom")
    h5ls(smallh5)
    h5read(smallh5, name = "Nhidden")
    smallsce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(logcounts = h5read(smallh5, name = "Y"))
    )
    #model_small <- newSlalomModel(smallsce)

})


test_that("SlalomModel trains", {
    gmtfile <- system.file("extdata", "reactome_subset.gmt", package = "slalom")
    genesets <- GSEABase::getGmt(gmtfile)
    data("mesc")
    model <- newSlalomModel(mesc, genesets, n_hidden = 5, min_genes = 10)
    model <- init(model, seed = 2222)
    expect_that(model, is_a("Rcpp_SlalomModel"))

    model <- train(model, nIterations = 5)

    library(rhdf5)
    h5file <- system.file("extdata", "outSim", "G_500N_20", "data.h5py",
                          package = "slalom")
    h5ls(h5file)
    h5read(h5file, "Nhidden")
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(logcounts = h5read(h5file, "Y"))
    )
    rownames(sce) <- unique(unlist(GSEABase::geneIds(genesets[1:20])))[1:500]
    model20 <- newSlalomModel(sce, genesets[1:20], n_hidden = 2, min_genes = 1)
    expect_error(init(model20, pi_prior = t(h5read(h5file, "Pi"))),
                      "n_hidden")
    model20 <- init(model20, pi_prior = t(h5read(h5file, "Pi")), n_hidden = 2)

    model20 <- train(model20, nIterations = 10)

    m1 <- newSlalomModel(sce, genesets[1:20], n_hidden = 2, min_genes = 1)
    m1 <- init(m1, pi_prior = t(h5read(h5file, "Pi")), n_hidden = 2)
    m1 <- slalom:::updateSlalom(m1)
    Zr <- m1$X_E1 %*% t(m1$W_E1)
    Zd <- m1$Z_E1 - Zr
    mean(abs(Zd))
    m1$Z_E1 <- m1$X_E1 %*% t(m1$W_E1)

    m2 <- newSlalomModel(sce, genesets[1:20], n_hidden = 2, min_genes = 1)
    m2 <- init(m2, pi_prior = t(h5read(h5file, "Pi")), n_hidden = 2)
    m2 <- slalom:::updateSlalom(m2)
    Zr <- m2$X_E1 %*% t(m2$W_E1)
    Zd <- m2$Z_E1 - Zr
    mean(abs(Zd))
    m2$Z_E1 <- m2$X_E1 %*% t(m2$W_E1)
    m2 <- slalom:::updateSlalom(m2)
    Zr <- m2$X_E1 %*% t(m2$W_E1)
    Zd <- m2$Z_E1 - Zr
    mean(abs(Zd))
    any(m2$epsilon_E1 < 0)
    sum(m2$epsilon_E1 < 0)
    any(is.nan(m2$X_E1))
    any(is.nan(m2$W_E1))
    m2$Z_E1 <- m2$X_E1 %*% t(m2$W_E1)
    m2 <- slalom:::updateSlalom(m2)
    any(m2$epsilon_E1 < 0)
    sum(m2$epsilon_E1 < 0)
    any(is.nan(m2$X_E1))
    any(is.nan(m2$W_E1))
    Zr <- m2$X_E1 %*% t(m2$W_E1)
    Zd <- m2$Z_E1 - Zr
    mean(abs(Zd))

    h5file2 <- system.file("extdata", "outSim", "20_cells_res_20it.h5py",
                           package = "slalom")
    h5ls(h5file2)
    h5read(h5file2, "Alpha")
    h5read(h5file2, "Pi_inferred")[1,,]
    dim(h5read(h5file2, "Pi_inferred")[20,,])
    dim(h5read(h5file2, "I"))
    h5read(h5file2, "I")[1:10, 1:7]

    m1 <- newSlalomModel(sce, genesets[1:20], n_hidden = 2, min_genes = 1)
    m1 <- init(m1, pi_prior = h5read(h5file2, "Pi_inferred")[1,,], n_hidden = 2)
    m1$X_E1 <- h5read(h5file2, "X")[1,,]
    m1$W_E1 <- h5read(h5file2, "W")[1,,]
    m1$alpha_E1 <- h5read(h5file2, "Alpha")[1,]
    m1 <- slalom:::updateSlalom(m1)
    expect_equal(m1$Pi_E1, h5read(h5file2, "Pi_inferred")[1,,])
    expect_equal(m1$X_E1, h5read(h5file2, "X")[2,,])
    View(m1$X_E1)
    expect_equal(m1$W_E1, h5read(h5file2, "W")[2,,])

    mem_change(m1 <- train(m1, nIterations = 800))
    plot(m1$W_E1, h5read(h5file2, "W")[21,,])

    m21 <- newSlalomModel(sce, genesets[1:20], n_hidden = 2, min_genes = 1)
    m21 <- init(m21, pi_prior = h5read(h5file2, "Pi_inferred")[21,,], n_hidden = 2)
    m21$X_E1 <- h5read(h5file2, "X")[21,,]
    m21$W_E1 <- h5read(h5file2, "W")[21,,]
    m21$alpha_E1 <- h5read(h5file2, "Alpha")[21,]
    m21 <- slalom:::updateSlalom(m21)
    expect_equal(m21$Pi_E1, h5read(h5file2, "Pi_inferred")[21,,])
    expect_equal(m21$W_E1, h5read(h5file2, "W")[20,,])
    expect_equal(m21$X_E1, h5read(h5file2, "X")[20,,])
    View(m21$X_E1)
    m21$Z_E1 <- m21$X_E1 %*% t(m21$W_E1)
    m21 <- slalom:::updateSlalom(m21)


    m21 <- train(m21, nIterations = 10)
    m21 <- train(m21, nIterations = 800)
    plot(m21$W_E1, h5read(h5file2, "W")[1,,])

    m2$updateW(0)
    any(is.nan(m2$W_E1))
    any(is.nan(m2$W_gamma0))
    m2$updateW(1)
    any(is.nan(m2$W_E1))
    any(is.nan(m2$W_gamma0))
    m2$updateW(2)
    any(is.nan(m2$W_E1))
    any(is.nan(m2$W_gamma0))
    for (i in 1:21)
        m2$updateW(i)
    View(m2$W_E1)




})


