# tests for SlalomModel methods

library(slalom)

test_that("newSlalomModel produces a valid object", {
    gmtfile <- system.file("extdata", "reactome_subset.gmt", package = "slalom")
    genesets <- GSEABase::getGmt(gmtfile)
    data("mesc")
    model <- newSlalomModel(mesc, genesets, n_hidden = 5, min_genes = 10)
    expect_that(model, is_a("Rcpp_SlalomModel"))

    lib_size <- rnbinom(ncol(mesc), mu = 50000, size = 1)
    design <- model.matrix(~lib_size)
    model <- newSlalomModel(mesc, genesets, n_hidden = 5, min_genes = 10,
                            design = design)
    expect_that(model, is_a("Rcpp_SlalomModel"))

    })


test_that("SlalomModel initialises", {
    gmtfile <- system.file("extdata", "reactome_subset.gmt", package = "slalom")
    genesets <- GSEABase::getGmt(gmtfile)
    data("mesc")
    model <- newSlalomModel(mesc, genesets, n_hidden = 5, min_genes = 10)
    model <- initSlalom(model)
    expect_that(model, is_a("Rcpp_SlalomModel"))

})


test_that("SlalomModel trains", {
    gmtfile <- system.file("extdata", "reactome_subset.gmt", package = "slalom")
    genesets <- GSEABase::getGmt(gmtfile)
    # for (i in seq_along(genesets)) {
    #     GSEABase::setName(genesets[[i]]) <- gsub("REACTOME_", "",
    #                                              GSEABase::setName(genesets[[i]]))
    #     GSEABase::setName(genesets[[i]]) <- strtrim(GSEABase::setName(genesets[[i]]), 30)
    # }
    ###########################################################################
    ## simulated data
    rdsfile <- system.file("extdata", "sim_N_20_v3.rds", package = "slalom")
    sim <- readRDS(rdsfile)
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(logcounts = sim[["init"]][["Y"]])
    )
    rownames(sce) <- unique(unlist(GSEABase::geneIds(genesets[1:20])))[1:500]

    ## test R and Py results with exactly the same initialisation
    m <- newSlalomModel(sce, genesets[1:23], n_hidden = 1, min_genes = 1)
    ## initialise this model
    m <- initSlalom(m, pi_prior = sim[["init"]][["Pi"]], n_hidden = 1)
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
    m <- trainSlalom(m, minIterations = 400, nIterations = 500, shuffle = FALSE,
                     pretrain = FALSE)
    expect_true(cor(c(m$alpha_E1[,1]),
                    c(sim[["final_iter"]][["Alpha"]])) > 0.9999)
    expect_true(cor(c(m$X_E1), c(sim[["final_iter"]][["X"]])) > 0.999999)
    expect_true(cor(c(m$W_E1), c(sim[["final_iter"]][["W"]])) > 0.999999)
    expect_true(cor(c(m$epsilon_E1[,1]),
                    c(sim[["final_iter"]][["Epsilon"]])) > 0.999999)
    expect_true(cor(c(m$W_gamma0), c(sim[["final_iter"]][["W_gamma0"]])) >
                    0.999999)

    ## test with pretraining and shuffling on
    m <- newSlalomModel(sce, genesets[1:23], n_hidden = 1, min_genes = 1)
    ## initialise this model
    m <- initSlalom(m, pi_prior = sim[["init"]][["Pi"]], n_hidden = 1, seed = 222)
    m$X_E1 <- sim[["init"]][["X"]]
    m$W_E1 <- sim[["init"]][["W"]]
    m$alpha_E1 <- sim[["init"]][["Alpha"]]
    m$epsilon_E1 <- sim[["init"]][["Epsilon"]]
    m$W_gamma0 <- sim[["init"]][["W_gamma0"]]
    m <- trainSlalom(m, minIterations = 400, nIterations = 500, shuffle = TRUE,
                     pretrain = TRUE, seed = 222)
    off_factors <- c(11, 12, 16, 18, 21)
    expect_true(cor(c(m$alpha_E1[,1]),
                    c(sim[["final_iter"]][["Alpha"]])) > 0.9)
    ## test corr of absolute values for all factors
    expect_true(cor(abs(c(m$X_E1)), abs(c(sim[["final_iter"]][["X"]]))) > 0.90)
    ## test corr of factors that are "on"
    expect_true(cor(c(m$X_E1[, 2:6]), c(sim[["final_iter"]][["X"]][, 2:6])) >
                    0.999)
    ## test corr of absolute values for all factors
    expect_true(cor(abs(c(m$W_E1)), abs(c(sim[["final_iter"]][["W"]]))) > 0.99)
    expect_true(cor(abs(c(m$W_E1[, -off_factors])),
                    abs(c(sim[["final_iter"]][["W"]][, -off_factors]))) > 0.99)
    ## test corr of factors that are "on"
    expect_true(cor(c(m$W_E1[, 2:5]), c(sim[["final_iter"]][["W"]][, 2:5])) >
                    0.999)
    expect_true(cor(c(m$epsilon_E1[,1]),
                    c(sim[["final_iter"]][["Epsilon"]])) > 0.9)
    expect_true(cor(c(m$W_gamma0), c(sim[["final_iter"]][["W_gamma0"]])) >
                    0.998)
    ## run this model to convergence
    m <- newSlalomModel(sce, genesets[1:23], n_hidden = 1, min_genes = 1)
    ## initialise this model
    m <- initSlalom(m, pi_prior = sim[["init"]][["Pi"]], n_hidden = 1, seed = 222)
    m$X_E1 <- sim[["init"]][["X"]]
    m$W_E1 <- sim[["init"]][["W"]]
    m$alpha_E1 <- sim[["init"]][["Alpha"]]
    m$epsilon_E1 <- sim[["init"]][["Epsilon"]]
    m$W_gamma0 <- sim[["init"]][["W_gamma0"]]
    m <- trainSlalom(m, minIterations = 400, nIterations = 3000, shuffle = TRUE,
                     pretrain = TRUE, seed = 222)
    expect_true(m$converged)

    # test without using same initialisation as Python code
    mm <- newSlalomModel(sce, genesets[1:23], n_hidden = 1, min_genes = 1)
    ## initialise this model
    mm <- initSlalom(mm, pi_prior = sim[["init"]][["Pi"]], n_hidden = 1, seed = 222)
    mm <- trainSlalom(mm, minIterations = 400, nIterations = 5000, shuffle = TRUE,
                     pretrain = TRUE, seed = 222)
    expect_true(mm$converged)
    off_factors <- c(11, 12, 16, 17, 18, 21)
    Zm <- m$X_E1 %*% t(m$W_E1)
    Zmm <- mm$X_E1 %*% t(mm$W_E1)
    expect_true(cor(c(Zm), c(Zmm)) > 0.9999)

    # genes_corr <- rep(NA, ncol(Zm))
    # for (i in 1:ncol(Zm))
    #     genes_corr[i] <- cor(Zm[,i], Zmm[,i])
    # summary(genes_corr)
    #
    # par(mfrow = c(1, 1), mar = c(3.1, 2.1, 2.1, 2.1))
    # plot(Zm, Zmm)
    # abline(0, 1, col = "firebrick")
    # Zm <- m$X_E1[, -off_factors] %*% t(m$W_E1[, -off_factors])
    # Zmm <- mm$X_E1[, -off_factors] %*% t(mm$W_E1[, -off_factors])
    # expect_true(cor(c(Zm), c(Zmm)) > 0.9999)
    # cor(c(Zm), c(Zmm))
    # par(mfrow = c(3, 7), mar = c(0.1, 0.1, 0.1, 0.1))
    # for (i in 1:20) {
    #     plot(Zm[i,], Zmm[i,])
    #     text(-1, 0, labels = i)
    #     abline(0, 1, col = "firebrick")
    # }
    # for (i in 1:20) {
    #     plot(Zm[i,], mm$Y[i,])
    #     text(-1, 0, labels = i)
    #     abline(0, 1, col = "firebrick")
    # }
    # for (i in 1:20) {
    #     plot(Zmm[i,], mm$Y[i,])
    #     text(-1, 0, labels = i)
    #     abline(0, 1, col = "firebrick")
    # }
    # par(mfrow = c(3, 7), mar = c(0.1, 0.1, 0.1, 0.1))
    # for (i in 1:21) {
    #     plot(m$W_E1[,i], mm$W_E1[,i])
    #     text(-1, 0, labels = i)
    #     abline(0, 1, col = "firebrick")
    # }
    # for (i in 1:21) {
    #     plot(m$W_gamma0[,i], mm$W_gamma0[,i])
    #     text(-1, 0, labels = i)
    #     abline(0, 1, col = "firebrick")
    # }
    # for (i in 1:21) {
    #     plot(m$X_E1[,i], mm$X_E1[,i])
    #     text(-1, 0, labels = i)
    #     abline(0, 1, col = "firebrick")
    # }
    # par(mfcol = c(1, 2))
    # plot(m$alpha_E1[,1], mm$alpha_E1[,1])
    # abline(0, 1, col = "firebrick")
    # plot(m$epsilon_E1[,1], mm$epsilon_E1[,1])
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
    model <- initSlalom(model, seed = 222)
    expect_that(model, is_a("Rcpp_SlalomModel"))

    model <- trainSlalom(model, nIterations = 2000)
    model$alpha_E1

    expect_that(model, is_a("Rcpp_SlalomModel"))
    expect_true(model$converged)

    # model <- newSlalomModel(mesc, genesets[1:100], n_hidden = 5, min_genes = 10,
    #                         prune_genes = FALSE)
    # model <- initSlalom(model, seed = 222)
    # expect_that(model, is_a("Rcpp_SlalomModel"))
    # model <- trainSlalom(model, nIterations = 2000)
    # model$alpha_E1
    #
    # expect_that(model, is_a("Rcpp_SlalomModel"))
    # expect_true(model$converged)

    # model <- newSlalomModel(mesc, genesets, n_hidden = 5, min_genes = 10)
    # model <- initSlalom(model)
    # expect_that(model, is_a("Rcpp_SlalomModel"))
    #
    # model <- trainSlalom(model, nIterations = 3000)
    #
    # expect_that(model, is_a("Rcpp_SlalomModel"))
    # expect_true(model$converged)

})


test_that("SlalomModel works with covariates", {
    # gmtfile <- system.file("extdata", "reactome_subset.gmt", package = "slalom")
    # genesets <- GSEABase::getGmt(gmtfile)
    # names(genesets) <- gsub("REACTOME_", "", names(genesets))
    # names(genesets) <- strtrim(names(genesets), 20)
    # ###########################################################################
    # ## simulated data
    # rdsfile <- system.file("extdata", "sim_N_20_v3.rds", package = "slalom")
    # sim <- readRDS(rdsfile)
    # sce <- SingleCellExperiment::SingleCellExperiment(
    #     assays = list(logcounts = sim[["init"]][["Y"]])
    # )
    # rownames(sce) <- unique(unlist(GSEABase::geneIds(genesets[1:20])))[1:500]
    #
    # ## access Python results if needed
    # library(rhdf5)
    # h5file <- file.path("./inst", "extdata", "outSim", "N_20_20it_cov2.h5py")
    # h5ls(h5file)
    # View(h5read(h5file, "Xcov"))
    # View(t(h5read(h5file, "X")[,,1]))
    # View(t(h5read(h5file, "I")))
    # View(t(h5read(h5file, "Pi_inferred")[,,1]))
    #
    # sce <- SingleCellExperiment::SingleCellExperiment(
    #     assays = list(logcounts = h5read(h5file, "Y"))
    # )
    # rownames(sce) <- unique(unlist(GSEABase::geneIds(genesets[1:20])))[1:500]
    #
    # ## test R and Py results with exactly the same initialisation
    # m <- newSlalomModel(sce, genesets[1:23], n_hidden = 1, min_genes = 1,
    #                     design = matrix(h5read(h5file, "Xcov"), ncol = 1))
    # ## initialise this model
    # m <- initSlalom(m, pi_prior = t(h5read(h5file, "Pi_inferred")[,,1]),
    #                 n_hidden = 1, design = matrix(h5read(h5file, "Xcov")))
    # m$X_E1 <- t(h5read(h5file, "X")[,,1])
    # m$W_E1 <- t(h5read(h5file, "W")[,,1])
    # m$alpha_E1 <- h5read(h5file, "Alpha")[,1]
    # m$epsilon_E1 <- h5read(h5file, "tau")[,1]
    # expect_that(m, is_a("Rcpp_SlalomModel"))
    #
    # m$updateW(0)
    # m$updateAlpha(0)
    # m$updateX(0)
    # cbind(h5read(h5file, "Alpha")[,2], m$alpha_E1[,1])
    # par(mfrow = c(2, 2), mar = c(2.1, 2.1, 2.1, 2.1))
    # plot(t(h5read(h5file, "W")[,,2])[,1], m$W_E1[,1], main = "W_E1")
    # abline(0, 1, col = "firebrick")
    # plot(t(h5read(h5file, "Pi_inferred")[,,2])[,1], m$W_gamma0[,1], main = "W_gamma0")
    # abline(0, 1, col = "firebrick")
    # plot(t(h5read(h5file, "X")[,,1])[,1], m$X_E1[,1], main = "X_E1")
    # abline(0, 1, col = "firebrick")
    #
    #
    #
    # par(mfrow = c(4, 6), mar = c(0.1, 0.1, 0.1, 0.1))
    # for (i in 1:23) {
    #     plot(t(h5read(h5file, "W")[,,1])[,i], m$W_E1[,i])
    #     text(-1, 0, labels = i)
    #     abline(0, 1, col = "firebrick")
    # }
    # par(mfrow = c(4, 6), mar = c(0.1, 0.1, 0.1, 0.1))
    # for (i in 1:23) {
    #     plot(t(h5read(h5file, "Pi_inferred")[,,1])[,i], m$W_gamma0[,i])
    #     text(-1, 0, labels = i)
    #     abline(0, 1, col = "firebrick")
    # }
    # par(mfrow = c(4, 6), mar = c(0.1, 0.1, 0.1, 0.1))
    # for (i in 1:23) {
    #     plot(t(h5read(h5file, "X")[,,1])[,i], m$X_E1[,i])
    #     text(-1, 0, labels = i)
    #     abline(0, 1, col = "firebrick")
    # }
    # par(mfcol = c(1, 2))
    # plot(h5read(h5file, "Alpha")[,1], m$alpha_E1[,1])
    # abline(0, 1, col = "firebrick")
    # plot(h5read(h5file, "tau")[,1], m$epsilon_E1[,1])
    # abline(0, 1, col = "firebrick")
    # par(mfcol = c(1, 1))
    # plot(t(h5read(h5file, "Y")), m$X_E1 %*% t(m$W_E1))
    # abline(0, 1, col = "firebrick")
    #
    # ## test precise results from one update
    # m <- updateSlalom(m)
    # expect_equal(m$alpha_E1[,1], h5read(h5file, "Alpha")[,2])
    # expect_equal(m$X_E1, t(h5read(h5file, "X")[,,2]))
    # expect_equal(m$W_E1, t(h5read(h5file, "W")[,,2]))
    # expect_equal(m$epsilon_E1[,1], h5read(h5file, "tau")[,2])
    # expect_equal(m$W_gamma0, t(h5read(h5file, "Pi_inferred")[,,2]))
    #
    # cbind(h5read(h5file, "Alpha")[,2], m$alpha_E1[,1])
    # par(mfrow = c(4, 6), mar = c(0.1, 0.1, 0.1, 0.1))
    # for (i in 1:23) {
    #     plot(t(h5read(h5file, "W")[,,2])[,i], m$W_E1[,i])
    #     text(-1, 0, labels = i)
    #     abline(0, 1, col = "firebrick")
    # }
    # par(mfrow = c(4, 6), mar = c(0.1, 0.1, 0.1, 0.1))
    # for (i in 1:23) {
    #     plot(t(h5read(h5file, "Pi_inferred")[,,2])[,i], m$W_gamma0[,i])
    #     text(-1, 0, labels = i)
    #     abline(0, 1, col = "firebrick")
    # }
    # par(mfrow = c(4, 6), mar = c(0.1, 0.1, 0.1, 0.1))
    # for (i in 1:23) {
    #     plot(t(h5read(h5file, "X")[,,2])[,i], m$X_E1[,i])
    #     text(-1, 0, labels = i)
    #     abline(0, 1, col = "firebrick")
    # }
    # par(mfcol = c(1, 2))
    # plot(h5read(h5file, "Alpha")[,2], m$alpha_E1[,1])
    # abline(0, 1, col = "firebrick")
    # plot(h5read(h5file, "tau")[,2], m$epsilon_E1[,1])
    # abline(0, 1, col = "firebrick")
    # par(mfcol = c(1, 1))
    # plot(t(h5read(h5file, "Y")), m$X_E1 %*% t(m$W_E1))
    # abline(0, 1, col = "firebrick")
    #
    #
    # ## train for 500 iterations and compare to Python results
    # m <- trainSlalom(m, minIterations = 400, nIterations = 500, shuffle = FALSE,
    #                  pretrain = FALSE)
    #
    # par(mfrow = c(4, 6), mar = c(0.1, 0.1, 0.1, 0.1))
    # for (i in 1:23) {
    #     plot(t(h5read(h5file, "W")[,,22])[,i], m$W_E1[,i])
    #     text(-1, 0, labels = i)
    #     abline(0, 1, col = "firebrick")
    # }
    # par(mfrow = c(4, 6), mar = c(0.1, 0.1, 0.1, 0.1))
    # for (i in 1:23) {
    #     plot(t(h5read(h5file, "Pi_inferred")[,,22])[,i], m$W_gamma0[,i])
    #     text(-1, 0, labels = i)
    #     abline(0, 1, col = "firebrick")
    # }
    # par(mfrow = c(4, 6), mar = c(0.1, 0.1, 0.1, 0.1))
    # for (i in 1:23) {
    #     plot(t(h5read(h5file, "X")[,,22])[,i], m$X_E1[,i])
    #     text(-1, 0, labels = i)
    #     abline(0, 1, col = "firebrick")
    # }
    # par(mfcol = c(1, 2))
    # cbind(h5read(h5file, "Alpha")[,22], m$alpha_E1[,1])
    # plot(h5read(h5file, "Alpha")[,22], m$alpha_E1[,1])
    # abline(0, 1, col = "firebrick")
    # plot(h5read(h5file, "tau")[,22], m$epsilon_E1[,1])
    # abline(0, 1, col = "firebrick")
    #
    #
    # expect_true(cor(c(m$alpha_E1[,1]),
    #                 c(sim[["final_iter"]][["Alpha"]])) > 0.9999)
    # expect_true(cor(c(m$X_E1), c(sim[["final_iter"]][["X"]])) > 0.999999)
    # expect_true(cor(c(m$W_E1), c(sim[["final_iter"]][["W"]])) > 0.999999)
    # expect_true(cor(c(m$epsilon_E1[,1]),
    #                 c(sim[["final_iter"]][["Epsilon"]])) > 0.999999)
    # expect_true(cor(c(m$W_gamma0), c(sim[["final_iter"]][["W_gamma0"]])) >
    #                 0.999999)
    #
    # ## test with pretraining and shuffling on
    # m <- newSlalomModel(sce, genesets[1:20], n_hidden = 1, min_genes = 1)
    # ## initialise this model
    # m <- initSlalom(m, pi_prior = sim[["init"]][["Pi"]], n_hidden = 1, seed = 222)
    # m$X_E1 <- sim[["init"]][["X"]]
    # m$W_E1 <- sim[["init"]][["W"]]
    # m$alpha_E1 <- sim[["init"]][["Alpha"]]
    # m$epsilon_E1 <- sim[["init"]][["Epsilon"]]
    # m$W_gamma0 <- sim[["init"]][["W_gamma0"]]
    # m <- trainSlalom(m, minIterations = 400, nIterations = 500, shuffle = TRUE,
    #                  pretrain = TRUE, seed = 222)
    # off_factors <- c(11, 12, 16, 18, 21)
    # expect_true(cor(c(m$alpha_E1[,1]),
    #                 c(sim[["final_iter"]][["Alpha"]])) > 0.95)
    # ## test corr of absolute values for all factors
    # expect_true(cor(abs(c(m$X_E1)), abs(c(sim[["final_iter"]][["X"]]))) > 0.90)
    # ## test corr of factors that are "on"
    # expect_true(cor(c(m$X_E1[, 2:6]), c(sim[["final_iter"]][["X"]][, 2:6])) >
    #                 0.999)
    # ## test corr of absolute values for all factors
    # expect_true(cor(abs(c(m$W_E1)), abs(c(sim[["final_iter"]][["W"]]))) > 0.99)
    # expect_true(cor(abs(c(m$W_E1[, -off_factors])),
    #                 abs(c(sim[["final_iter"]][["W"]][, -off_factors]))) > 0.99)
    # ## test corr of factors that are "on"
    # expect_true(cor(c(m$W_E1[, 2:5]), c(sim[["final_iter"]][["W"]][, 2:5])) >
    #                 0.999)
    # expect_true(cor(c(m$epsilon_E1[,1]),
    #                 c(sim[["final_iter"]][["Epsilon"]])) > 0.93)
    # expect_true(cor(c(m$W_gamma0), c(sim[["final_iter"]][["W_gamma0"]])) >
    #                 0.998)
    # ## run this model to convergence
    # m <- newSlalomModel(sce, genesets[1:20], n_hidden = 1, min_genes = 1)
    # ## initialise this model
    # m <- initSlalom(m, pi_prior = sim[["init"]][["Pi"]], n_hidden = 1, seed = 222)
    # m$X_E1 <- sim[["init"]][["X"]]
    # m$W_E1 <- sim[["init"]][["W"]]
    # m$alpha_E1 <- sim[["init"]][["Alpha"]]
    # m$epsilon_E1 <- sim[["init"]][["Epsilon"]]
    # m$W_gamma0 <- sim[["init"]][["W_gamma0"]]
    # m <- trainSlalom(m, minIterations = 400, nIterations = 3000, shuffle = TRUE,
    #                  pretrain = TRUE, seed = 222)
    # expect_true(m$converged)
    #
    # # test without using same initialisation as Python code
    # mm <- newSlalomModel(sce, genesets[1:20], n_hidden = 1, min_genes = 1)
    # ## initialise this model
    # mm <- initSlalom(mm, pi_prior = sim[["init"]][["Pi"]], n_hidden = 1, seed = 222)
    # mm <- trainSlalom(mm, minIterations = 400, nIterations = 3000, shuffle = TRUE,
    #                   pretrain = TRUE, seed = 222)
    # expect_true(mm$converged)
    # off_factors <- c(11, 12, 16, 17, 18, 21)
    # Zm <- m$X_E1 %*% t(m$W_E1)
    # Zmm <- mm$X_E1 %*% t(mm$W_E1)
    # expect_true(cor(c(Zm), c(Zmm)) > 0.9999)
    #
    # genes_corr <- rep(NA, ncol(Zm))
    # for (i in 1:ncol(Zm))
    #     genes_corr[i] <- cor(Zm[,i], Zmm[,i])
    # summary(genes_corr)
    #
    # par(mfrow = c(1, 1), mar = c(3.1, 2.1, 2.1, 2.1))
    # plot(Zm, Zmm)
    # abline(0, 1, col = "firebrick")
    # Zm <- m$X_E1[, -off_factors] %*% t(m$W_E1[, -off_factors])
    # Zmm <- mm$X_E1[, -off_factors] %*% t(mm$W_E1[, -off_factors])
    # expect_true(cor(c(Zm), c(Zmm)) > 0.9999)
    # cor(c(Zm), c(Zmm))
    # par(mfrow = c(3, 7), mar = c(0.1, 0.1, 0.1, 0.1))
    # for (i in 1:20) {
    #     plot(Zm[i,], Zmm[i,])
    #     text(-1, 0, labels = i)
    #     abline(0, 1, col = "firebrick")
    # }
    # for (i in 1:20) {
    #     plot(Zm[i,], mm$Y[i,])
    #     text(-1, 0, labels = i)
    #     abline(0, 1, col = "firebrick")
    # }
    # for (i in 1:20) {
    #     plot(Zmm[i,], mm$Y[i,])
    #     text(-1, 0, labels = i)
    #     abline(0, 1, col = "firebrick")
    # }
    # par(mfrow = c(3, 7), mar = c(0.1, 0.1, 0.1, 0.1))
    # for (i in 1:21) {
    #     plot(m$W_E1[,i], mm$W_E1[,i])
    #     text(-1, 0, labels = i)
    #     abline(0, 1, col = "firebrick")
    # }
    # for (i in 1:21) {
    #     plot(m$W_gamma0[,i], mm$W_gamma0[,i])
    #     text(-1, 0, labels = i)
    #     abline(0, 1, col = "firebrick")
    # }
    # for (i in 1:21) {
    #     plot(m$X_E1[,i], mm$X_E1[,i])
    #     text(-1, 0, labels = i)
    #     abline(0, 1, col = "firebrick")
    # }
    # par(mfcol = c(1, 2))
    # plot(m$alpha_E1[,1], mm$alpha_E1[,1])
    # abline(0, 1, col = "firebrick")
    # plot(m$epsilon_E1[,1], mm$epsilon_E1[,1])
    # abline(0, 1, col = "firebrick")
    #



})


