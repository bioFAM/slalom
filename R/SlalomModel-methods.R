# Methods for the (Rcpp)SlalomModel class

################################################################################
### New SlalomModel class object

#' Create a new SlalomModel object.
#'
#' Slalom fits relatively complicated hierarchical Bayesian factor analysis
#' models with data and results stored in a \code{"SlalomModel"} object. This
#' function builds a new \code{"SlalomModel"} object from minimal inputs.
#'
#' @param object \code{"SingleCellExperiment"} object N x G expression data matrix (cells x genes)
#' @param genesets a \code{"GeneSetCollection"} object containing annotated
#' gene sets
#' @param n_hidden number of hidden factors to fit in the model (2-5 recommended)
#' @param prune_genes logical, should genes that are not annotated to any gene
#' sets be filtered out?
#' @param min_genes scalar, minimum number of genes required in order to retain
#' a gene set for analysis
#' @param design numeric design matrix providing values for covariates to fit in
#' the model (rows represent cells)
#' @param anno_fpr numeric(1), false positive rate (FPR) for assigning genes to
#' factors (pathways); default is 0.01
#' @param anno_fnr numeric(1), false negative rate (FNR) for assigning genes to
#' factors (pathways); default is 0.001
#' @param assay_name character(1), the name of the \code{assay} of the
#' \code{object} to use as expression values. Default is \code{logcounts},
#' assumed to be normalised log2-counts-per-million values or equivalent.
#' @param verbose logical(1), should information about what's going be printed
#' to screen?
#'
#' @return a new Rcpp_SlalomModel object
#'
#' @details
#' This function builds and returns the object, checking for validity, which
#' includes checking that the input data is of consistent dimensions.
#'
#' @importFrom GSEABase geneIds
#' @importFrom GSEABase getGmt
#' @import Rcpp
#' @import RcppArmadillo
#' @importFrom methods new
#' @importFrom methods validObject
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @export
#'
#' @examples
#' gmtfile <- system.file("extdata", "reactome_subset.gmt", package = "slalom")
#' genesets <- GSEABase::getGmt(gmtfile)
#' data("mesc")
#' model <- newSlalomModel(mesc, genesets, n_hidden = 5, min_genes = 10)
#'
#' exprsfile <- system.file("extdata", "mesc.csv", package = "slalom")
#' mesc_mat <- as.matrix(read.csv(exprsfile))
#' sce <- SingleCellExperiment::SingleCellExperiment(assays = list(logcounts = mesc_mat))
#' # model2 <- newSlalomModel(mesc_mat, genesets, n_hidden = 5, min_genes = 10)
#'
newSlalomModel <- function(
    object, genesets, n_hidden = 5, prune_genes = TRUE, min_genes = 15,
    design = NULL, anno_fpr = 0.01, anno_fnr = 0.001, assay_name = "logcounts",
    verbose = TRUE) {

    if (!methods::is(object, "SingleCellExperiment"))
        stop("object must be a SingleCellExperiment")
    Y <- t(SummarizedExperiment::assay(object, assay_name))
    if (is.null(Y))
        stop("object must contain 'logcounts' in assays")
    if (is.null(rownames(object)))
        stop("rownames(object) is NULL: expecting gene identifiers as rownames")
    colnames(Y) <- rownames(object)
    if (is.null(colnames(object)))
        rownames(Y) <- as.character(1:ncol(object))
    else
        rownames(Y) <- colnames(object)
    ## convert GeneSetCollection into I matrix of indicators
    I <- lapply(genesets, function(x) {colnames(Y) %in% GSEABase::geneIds(x)})
    names(I) <- names(genesets)
    I <- as.matrix(as.data.frame(I))
    rownames(I) <- colnames(Y)
    ## filter genesets on min_genes
    n_sets_pass <- colSums(I) >= min_genes
    if (sum(n_sets_pass) > 0L) {
        I <- I[, n_sets_pass]
        if (verbose)
            cat(sum(n_sets_pass), "annotated factors retained; ",
                length(genesets) - sum(n_sets_pass), "annotated factors dropped.\n")
    } else
        stop("No gene sets found containing more than min_genes genes.")
    ## prune genes if desired
    if (prune_genes) {
        keep_gene <- (rowSums(I) > 0L)
        I <- I[keep_gene, ]
    }
    retained_genes <- rownames(I)
    if (verbose)
        cat(length(retained_genes), " genes retained for analysis.\n")
    ## add hidden factors
    if (n_hidden > 0L) {
        hidden_factors <- matrix(TRUE, nrow = nrow(I), ncol = n_hidden)
        colnames(hidden_factors) <- paste0("hidden", sprintf("%02d", 1:n_hidden))
        I <- cbind(hidden_factors, I)
    }

    ## create new SlalomModel object for output
    ## set some attributes that we need frequently for the updates, inculuding
    ## number and idx of hidden (unannotated) terms
    pi_init <- I * 1L
    if (!is.null(design)) {
        if (is.null(colnames(design)))
            colnames(design) <- paste0("known", sprintf("%02d", 1:ncol(design)))
        tmpmat <- matrix(1, nrow = nrow(pi_init), ncol = ncol(design))
        colnames(tmpmat) <- colnames(design)
        pi_init <- cbind(tmpmat, pi_init)
    }
    pi_init[pi_init > 0.5] <- 1 - anno_fpr
    pi_init[pi_init < 0.5] <- anno_fnr
    out <- .newSlalom(Y[, retained_genes], pi_init, n_hidden, design)
    out$termNames <- c(colnames(design), colnames(I))
    out$cellNames <- rownames(Y)
    out$geneNames <- retained_genes
    out$pretrain_order <- seq_len(out$K) - 1
    ## Check validity of object
    #validObject(out)
    out
}


.newSlalom <- function(Y, pi_init, n_hidden, design = NULL) {
    n_known <- ifelse(is.null(design), 0, ncol(design))
    if (is.null(design))
        x_init <- matrix(0, nrow = nrow(Y), ncol = ncol(pi_init))
    else
        x_init <- cbind(
            design,
            matrix(0, nrow = nrow(Y), ncol = ncol(pi_init) - ncol(design)))
    Rcpp::loadModule("SlalomModel")
    out <- new(
        "Rcpp_SlalomModel",
        Y_init = Y,
        pi_init = pi_init,
        X_init = x_init,
        W_init = matrix(0, nrow = ncol(Y), ncol(pi_init)),
        prior_alpha = c(1e-03, 1e-03),
        prior_epsilon = c(1e-03, 1e-03)
    )
    ## add more data to the object
    out$Z_E1 <- out$X_E1 %*% t(out$W_E1 * pi_init)
    out$iUnannotatedDense <- seq(from = n_known + 1, to = n_known + n_hidden)
    out$nAnnotated <- ncol(pi_init) - n_hidden - n_known
    out$nHidden <- n_hidden
    out$nKnown <- n_known
    ## add design matrix as known factors if supplied
    if (!is.null(design)) {
        if (nrow(design) != out$N)
            stop("Number of rows of design does not match number of cells.")
        else {
            out$Known <- design
        }
    }
    out
}


################################################################################
### Initialize a SlalomModel class object

#' Initialize a SlalomModel object
#'
#' Initialize a SlalomModel with sensible starting values for parameters before
#' training the model.
#'
#' @param object a \code{Rcpp_SlalomModel} object
#' @param alpha_priors numeric(2) giving alpha and beta hyperparameters for a
#' gamma prior distribution for alpha parameters (precision of factor weights)
#' @param epsilon_priors numeric(2) giving alpha and beta hyperparameters for a
#' gamma prior distribution for noise precision parameters
#' @param noise_model character(1) defining noise model, defaults to "gauss" for
#' Gaussian noise model
#' @param seed integer(1) value supplying a random seed to make results
#' reproducible (default is \code{NULL})
#' @param pi_prior numeric matrix (genes x factors) giving prior probability of
#' a gene being active for a factor
#' @param n_hidden integer(1), number of hidden factors in model. Required if
#' \code{pi_prior} is not \code{NULL}, ignored otherwise.
#' @param design matrix of known factors (covariates) to fit in the
#' model. Optional if \code{pi_prior} is not \code{NULL}, ignored otherwise.
#' @param verbose logical(1), should messages be printed about what the function
#' is doing? Default is \code{TRUE}.
#' @param save_init logical(1), save the initial X values (factor states for
#' each cell) in the object? Default is \code{FALSE} as this is potentially a
#' large N (number of cell) by K (number of factors) matrix.
#'
#' @details It is strongly recommended to use \code{\link{newSlalomModel}} to
#' create the \code{\link{SlalomModel}} object prior to applying
#' \code{initSlalom}.
#'
#' @return an `Rcpp_SlalomModel` object
#'
#' @name initSlalom
#' @rdname initSlalom
#'
#' @author Davis McCarthy
#' @importFrom rsvd rpca
#' @importFrom stats runif
#' @importFrom stats rnorm
#' @importFrom stats prcomp
#' @export
#' @examples
#' gmtfile <- system.file("extdata", "reactome_subset.gmt", package = "slalom")
#' genesets <- GSEABase::getGmt(gmtfile)
#' data("mesc")
#' model <- newSlalomModel(mesc, genesets, n_hidden = 5, min_genes = 10)
#' model <- initSlalom(model)
initSlalom <- function(
    object, alpha_priors = NULL, epsilon_priors = NULL,  noise_model = "gauss",
    seed = NULL, pi_prior = NULL, n_hidden = NULL, design = NULL,
    verbose = FALSE, save_init = FALSE) {
    if (!methods::is(object, "Rcpp_SlalomModel"))
        stop("object must be of class Rcpp_SlalomModel")
    ## if Pi priors are supplied, need to redefine object
    if (!is.null(pi_prior)) {
        if (is.null(n_hidden))
            stop("If pi_prior is supplied, n_hidden and n_known must be provided.
                First n_known cols of pi_prior are taken as known covariates and next n_hidden cols are taken as hidden factors.")
        tnames <- object$termNames
        cnames <- object$cellNames
        gnames <- object$geneNames
        object <- .newSlalom(object$Y, pi_prior, n_hidden, design)
        object$termNames <- tnames
        object$cellNames <- cnames
        object$geneNames <- gnames
        object$pretrain_order <- seq_len(object$K) - 1
    }
    ## define noise model
    noise_model <- match.arg(noise_model, c("gauss")) ## future: support "hurdle", "poisson"
    object$noiseModel <- noise_model
    ## define priors for alpha and epsilon
    if (!is.null(alpha_priors)) {
        object$alpha_pa <- alpha_priors[0]
        object$alpha_pb <- alpha_priors[1]
        object$alpha_a <- rep(object$alpha_pa, object$K)
        object$alpha_b <- rep(object$alpha_pb, object$K)
        object$alpha_E1 <- object$alpha_b / object$alpha_a
        object$alpha_lnE1 <- digamma(object$alpha_a) - log(object$alpha_b)
    }
    if (!is.null(epsilon_priors)) {
        object$epsilon_pa <- epsilon_priors[0]
        object$epsilon_pb <- epsilon_priors[1]
        object$epsilon_a <- rep(object$epsilon_pa, object$G)
        object$epsilon_b <- rep(object$epsilon_pb, object$G)
        object$epsilon_E1 <- object$epsilon_b / object$epsilon_a
        object$epsilon_lnE1 <- digamma(object$epsilon_a) - log(object$epsilon_b)
    }

    ## pi is likelihood of link for genes x factors, defined in object from
    ## annotated gene sets with newSlalomModel
    ## define the number of genes that are "on" for each factor
    object$nOn <- colSums(object$Pi_E1 > 0.5)

    ## define initial expected values for indicator variables Z with observed
    ## annotations
    object$Z_E1 <- object$Y

    ## bits and bobs
    object$onF <- 1 / object$nScale # object@nScale
    # object@isExpressed <- (object@Z[["E1"]] > 0) * 1.0
    # object@numExpressed <- colSums(object@Z[["E1"]] > 0)

    ## initialize parameters using random PCA method
    object <- .initializePCARand(object, seed = seed, save_init = save_init)
    ## check validity of object and return
    validObject(object)
    object
}


.initializePCARand <- function(object, seed = NULL, save_init = FALSE) {
    if (!is.null(seed))
        set.seed(seed)
    # Zstd <- object$Z[["E1"]]
    ## initialize some objects
    Ystd <- object$Y
    unif_vars <- stats::runif(nrow(object$Pi_E1) * ncol(object$Pi_E1))
    Z_init <- object$Pi_E1
    Z_init[Z_init < 0.1] <- 0.1
    Z_init[Z_init > 0.9] <- 0.9
    Ion <- matrix(unif_vars, nrow = nrow(object$Pi_E1)) < Z_init
    object$X_E1 <- matrix(nrow = object$N, ncol = object$K)
    object$X_diagSigmaS <- matrix(1.0 / 2, nrow = object$N, ncol = object$K)
    object$W_E1 <- matrix(nrow = object$G, ncol = object$K)
    ## initialize gamma component of weights
    object$W_gamma0 <- object$Pi_E1
    object$W_gamma0[object$W_gamma0 <= 0.1] <- 0.1
    object$W_gamma0[object$W_gamma0 >= 0.9] <- 0.9
    object$W_gamma1 <- 1.0 - object$W_gamma0
    ## initialise known factors
    if (object$nKnown > 0L) {
        for (k in seq_len(object$nKnown)) {
            object$W_E1[, k] <- sqrt(1.0 / object$K) * stats::rnorm(object$G)
            object$X_diagSigmaS[, k] <- 1.0 / 2
        }
        object$X_E1[, seq_len(object$nKnown)] <- object$Known
    }
    ## initialise hidden (latent) factors
    if (object$nHidden > 0L) {
        for (iL in object$iUnannotatedDense) {
            object$X_E1[, iL] <- stats::rnorm(object$N)
            object$W_E1[, iL] <- sqrt(1.0 / object$K) * stats::rnorm(object$G)
        }
    }
    ## initialise annotated gene set factors
    for (k in seq_len(object$nAnnotated)) {
        k <- k + object$nKnown + object$nHidden
        object$W_E1[, k] <- sqrt(1.0 / object$K) * stats::rnorm(object$G)
        object$X_diagSigmaS[, k] <- 1.0 / 2
        if (sum(Ion[,k]) > 5) {
            if (object$N < 500)
                pca <- stats::prcomp(Ystd[, Ion[, k]], rank. = 1, retx = TRUE)
            else
                pca <- rsvd::rpca(Ystd[, Ion[, k]], k = 1, retx = TRUE)
            object$X_E1[, k] <- scale(pca$x[, 1])
        } else {
            object$X_E1[, k] <- stats::rnorm(object$N)
        }
    }
    ## save initial X values if desired
    if (save_init) {
        object$X_init <- object$X_E1
    }
    object
}



################################################################################
### Train a SlalomModel class object

#' Train a SlalomModel object
#'
#' Train a SlalomModel to infer model parameters.
#'
#' @param object a \code{Rcpp_SlalomModel} object
#' @param nIterations integer(1) maximum number of iterations to use in training
#' the model (default: 5000)
#' @param minIterations integer(1) minimum number of iterations to perform.
#' @param tolerance numeric(1) tolerance to allow between iterations
#' (default 1e-08)
#' @param forceIterations logical(1) should the model be forced to update
#' \code{nIteration} times?
#' @param shuffle logical(1) should the order in which factors are updated be
#' shuffled between iterations? Shuffling generally helps speed up convergence
#' so is recommended and defaults is \code{TRUE}
#' @param pretrain logical(1), should the model be "pre-trained" to achieve
#' faster convergence and obtain an initial update order? Recommended; default
#' is \code{TRUE}
#' @param verbose logical(1), should messages be printed about what the function
#' is doing? Default is \code{TRUE}.
#' @param seed integer(1) value supplying a random seed to make results
#' reproducible (default is \code{NULL})
#' @param drop_factors logical(1), should factors be dropped from the model if
#' the model determines them not to be relevant? Default is \code{TRUE}.
#'
#' @details Train the model using variational Bayes methods to infer parameters.
#'
#' @return an `Rcpp_SlalomModel` object
#'
#' @importFrom stats cor
#' @docType methods
#' @aliases train train,Rcpp_SlalomModel-method
#'
#' @author Davis McCarthy
#' @export
#' @examples
#' gmtfile <- system.file("extdata", "reactome_subset.gmt", package = "slalom")
#' genesets <- GSEABase::getGmt(gmtfile)
#' data("mesc")
#' model <- newSlalomModel(mesc, genesets, n_hidden = 5, min_genes = 10)
#' model <- initSlalom(model)
#' model <- trainSlalom(model, nIterations = 10)
trainSlalom <- function(
    object, nIterations = 5000, minIterations = 700, tolerance = 1e-08,
    forceIterations = FALSE, shuffle = TRUE, pretrain = TRUE, verbose = TRUE,
    seed = NULL, drop_factors = TRUE) {
    ## define training parameters in object
    object$dropFactors <- drop_factors  ## drop factors if deemed to be off
    object$tolerance <- tolerance
    object$nIterations <- nIterations
    object$minIterations <- minIterations
    object$forceIterations <- forceIterations
    object$shuffle <- shuffle ## shuffle order in which factors are updated each
    ## for each training iteration
    if (pretrain) {
        if (verbose)
            cat("pre-training model for faster convergence\n")
        pt_out <- .preTrain(object, seed = seed)
        object$pretrain_order <- pt_out - 1
    } else {
        object$pretrain_order <- seq_len(object$K) - 1
    }
    ## train model
    object$train()
    ## return object
    object
}


.preTrain <- function(object, seed = NULL, n_fix = NULL) {
    ## n_fix denotes the number of terms which should be fixed and updated first;
    ## defaults to NULL, which results in the nubmer of unannotated factors
    ## being updated first
    if (is.null(n_fix))
        n_fix <- object$nKnown + object$nHidden
    pi0 <- object$Pi_E1
    pi0[pi0 > 0.8] <- 1 - 1e-100
    pi0[pi0 < 0.2] <- 1e-100
    ## fit pca
    if (!is.null(seed))
        set.seed(seed)
    if (object$N < 500)
        pca <- stats::prcomp(object$Y, rank. = 1, retx = TRUE)
    else
        pca <- rsvd::rpca(object$Y, k = 1, retx = TRUE)
    ## sort by correlation to PC1
    mpc <- abs(stats::cor(object$X_E1, pca$x))[-c(1:n_fix)]
    Ipi <- order(mpc, decreasing = TRUE)
    IpiRev <- rev(Ipi)
    ## organise for fitting fwd and rev models
    k_range <- 1:object$K
    k_range[-c(1:n_fix)] <- Ipi + n_fix
    k_range_rev <- 1:object$K
    k_range_rev[-c(1:n_fix)] <- IpiRev + n_fix
    n_hidden <- object$nHidden
    if (object$nKnown > 0)
        design <- object$Known
    else
        design <- NULL
    Y <- object$Z_E1
    rm(object)
    ## run model for 50 iterations
    pi <- pi0[, k_range]
    obj_fwd <- .newSlalom(Y, pi, n_hidden, design = design)
    obj_fwd$shuffle <- TRUE
    obj_fwd$nScale <- 30
    obj_fwd <- initSlalom(obj_fwd)
    obj_fwd <- trainSlalom(obj_fwd, pretrain = FALSE, shuffle = FALSE,
                           nIterations = 50, minIterations = 50,
                           verbose = FALSE, drop_factors = FALSE)
    alpha_fwd <- obj_fwd$alpha_E1
    rm(obj_fwd)
    ## run reverse model for 50 iterations
    pi <- pi0[, k_range_rev]
    obj_rev <- .newSlalom(Y, pi, n_hidden, design = design)
    obj_rev$shuffle <- TRUE
    obj_rev$nScale <- 30
    obj_rev <- initSlalom(obj_rev)
    obj_rev <- trainSlalom(obj_rev, pretrain = FALSE, shuffle = FALSE,
                           nIterations = 50, minIterations = 50,
                           verbose = FALSE, drop_factors = FALSE)
    alpha_rev <- obj_rev$alpha_E1
    rm(obj_rev)
    ## get update ordering
    IpiK <- order(
        -(0.5 * (1.0 / alpha_rev[order(k_range_rev)][-c(1:n_fix)])) +
            0.5 * (1.0 / alpha_fwd[order(k_range)][-c(1:n_fix)])
    ) ## check this!!
    i_label <- c(1:n_fix, IpiK + n_fix)
    i_label
}



################################################################################
### Update a SlalomModel class object

#' Update a SlalomModel object
#'
#' Do one variational update of a SlalomModel to infer model parameters.
#'
#' @param object a \code{Rcpp_SlalomModel} object
#'
#' @details Update the model with one iteration using variational Bayes methods
#' to infer parameters.
#'
#' @name updateSlalom
#' @rdname updateSlalom
#' @importFrom methods is
#'
#' @return an `Rcpp_SlalomModel` object
#'
#' @author Davis McCarthy
#' @export
#' @examples
#' gmtfile <- system.file("extdata", "reactome_subset.gmt", package = "slalom")
#' genesets <- GSEABase::getGmt(gmtfile)
#' data("mesc")
#' model <- newSlalomModel(mesc, genesets, n_hidden = 5, min_genes = 10)
#' model <- initSlalom(model)
#' model <- updateSlalom(model)
updateSlalom <- function(object) {
    if (!methods::is(object, "Rcpp_SlalomModel"))
        stop("object must be of class Rcpp_SlalomModel")
    ## train model
    object$update()
    ## return object
    object
}


################################################################################
### Define validity check for SlalomModel class object

setValidity("Rcpp_SlalomModel", function(object) {
    msg <- NULL
    valid <- TRUE

    if ( is.null(object@Y) ) {
        valid <- FALSE
        msg <- c(msg,
                 "object must contain an expression matrix at object@Y")
    }
    if ( any(is.na(object@Y)) ) {
        valid <- FALSE
        msg <- c(msg,
                 "Expression matrix Y cannot contain NA values")
    }
    if ( any(is.na(object@Pi_E1)) ) {
        valid <- FALSE
        msg <- c(msg,
                 "Matrix Pi of observed gene set annotations cannot contain NA values")
    }
    if ( (ncol(object@Y) != object@G)) {
        valid <- FALSE
        msg <- c(msg,
                 "Number of genes in expression matrix Y must match object@G")
    }
    if ( (nrow(object@Y) != object@N)) {
        valid <- FALSE
        msg <- c(msg,
                 "Number of cells in expression matrix Y must match object@N")
    }
    if ( (ncol(object@Pi_E1) != object@K)) {
        valid <- FALSE
        msg <- c(msg,
                 "Number of factors in matrix pi must match object@K")
    }
    if ( (nrow(object@Pi_E1) != object@G)) {
        valid <- FALSE
        msg <- c(msg,
                 "Number of factors in matrix pi must match object@G")
    }
    if ( (length(object@termNames) != object@K)) {
        valid <- FALSE
        msg <- c(msg,
                 "Number of factors termNames must match object@K")
    }
    if ( (length(object@geneNames) != object@G)) {
        valid <- FALSE
        msg <- c(msg,
                 "Number of geneNames must match object@G")
    }
    if ( (length(object@cellNames) != object@N)) {
        valid <- FALSE
        msg <- c(msg,
                 "Number of cellNames must match object@N")
    }


    if (valid) TRUE else msg
})

