# Methods for the (Rcpp)SlalomModel class

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

#' @return a new SlalomModel object
#'
#' @details
#' This function builds and returns the object, checking for validity, which
#' includes checking that the input data is of consistent dimensions.
#'
#' @importFrom GSEABase geneIds
#' @importFrom GSEABase getGmt
#' @import scater
#' @import Rcpp
#' @import RcppArmadillo
#' @importFrom Biobase exprs
#' @importFrom methods new
#' @importFrom methods validObject
#' @importFrom Rcpp evalCpp
#' @useDynLib slalom
#' @aliases Rcpp_SlalomModel
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
newSlalomModel <- function(object, genesets, n_hidden = 5, prune_genes = TRUE,
                           min_genes = 15, design = NULL) {

    if (!is(object, "SingleCellExperiment"))
        stop("object must be a SingleCellExperiment")
    Y <- t(object@assays$data$logcounts)
    if (is.null(rownames(object)))
        stop("rownames(object) is NULL: expecting gene identifiers as rownames")
    if (is.null(colnames(Y)))
        colnames(Y) <- rownames(object)
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
    ## set some attributes that we need frequently for the updates, inculuding
    ## number and idx of hidden (unannotated) terms
    out <- .newSlalom(Y[, retained_genes], I * 1L, n_hidden, design)

    ## Check validity of object
    #validObject(out)
    out
}


.newSlalom <- function(Y, I, n_hidden, design = NULL) {
    n_known <- ifelse(is.null(design), 0, ncol(design))
    slalom_module <- Rcpp::Module("SlalomModel", PACKAGE = "slalom")
    out <- new(slalom_module$SlalomModel,
               Y_init = Y,
               pi_init = I,
               X_init = matrix(0, nrow = nrow(Y), ncol(I)),
               W_init = matrix(0, nrow = ncol(Y), ncol(I)),
               prior_alpha = c(1e-03, 1e-03),
               prior_epsilon = c(1e-03, 1e-03)
    )
    ## add more data to the object
    out$Z_E1 <- out$X_E1 %*% t(out$W_E1 * (1L * I))
    out$iUnannotatedDense <- seq(from = (ncol(I) - n_hidden + 1), to = ncol(I))
    out$nAnnotated <- ncol(I) - n_hidden
    out$nHidden <- n_hidden
    out$nKnown <- n_known
    ## add design matrix as known factors if supplied
    if (!is.null(design)) {
        if (nrow(design) != out$N)
            stop("Number of rows of design does not match number of cells.")
        else {
            out$K <- n_known + ncol(I)
            out$Known <- design
        }
    }
    out
}


# out <- new("SlalomModel",
#     K = n_known + ncol(I),
#     N = nrow(Y),
#     G = nrow(I),
#     nAnnotated = ncol(I) - n_hidden,
#     nHidden = n_hidden,
#     nKnown = n_known,
#     nScale = 100,
#     Y = Y[, retained_genes],
#     pi = I * 1L,
#     termNames = colnames(I),
#     geneNames = retained_genes,
#     cellNames = scater::cellNames(object),
#     iUnannotatedDense = seq(from = (ncol(I) - n_hidden + 1), to = ncol(I))
# )


################################################################################
### Initialize a SlalomModel class object

#' #' Initialize a SlalomModel object
#' #'
#' #' Initialize a SlalomModel with sensible starting values for parameters before
#' #' training the model.
#' #'
#' #' @param object a \code{Rcpp_SlalomModel} object
#' #' @param alpha_priors numeric(2) giving alpha and beta hyperparameters for a
#' #' gamma prior distribution for alpha parameters (precision of factor weights)
#' #' @param epsilon_priors numeric(2) giving alpha and beta hyperparameters for a
#' #' gamma prior distribution for noise precision parameters
#' #' @param noise_model character(1) defining noise model, defaults to "gauss" for
#' #' Gaussian noise model
#' #' @param seed integer(1) value supplying a random seed to make results
#' #' reproducible
#' #' @param ... generic arguments passed to \code{Rcpp_SlalomModel} method
#' #'
#' #' @details It is strongly recommended to use \code{\link{newSlalomModel}} to
#' #' create the \code{\link{SlalomModel}} object prior to applying
#' #' \code{initialize}.
#' #'
#' #' @docType methods
#' #' @name init
#' #' @rdname init
#' #' @aliases init init,Rcpp_SlalomModel-method
#' #'
#' #' @author Davis McCarthy
#' #' @importFrom rsvd rpca
#' #' @importFrom stats runif
#' #' @importFrom stats rnorm
#' #' @export
#' #' @examples
#' #' gmtfile <- system.file("extdata", "reactome_subset.gmt", package = "slalom")
#' #' genesets <- GSEABase::getGmt(gmtfile)
#' #' data("mesc")
#' #' model <- newSlalomModel(mesc, genesets, n_hidden = 5, min_genes = 10)
#' #' model <- init(model)
#' ## setGeneric("init", function(object, ...) standardGeneric("init"))
#'
#' #' @name init
#' #' @rdname init
#' #' @docType methods
#' #' #@export
#' ## setMethod("init", "Rcpp_SlalomModel", function(


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
#' reproducible
#' @param pi_prior numeric matrix (genes x factors) giving prior probability of
#' a gene being active for a factor
#'
#' @details It is strongly recommended to use \code{\link{newSlalomModel}} to
#' create the \code{\link{SlalomModel}} object prior to applying
#' \code{initialize}.
#'
#' @name init
#' @rdname init
#'
#' @author Davis McCarthy
#' @importFrom rsvd rpca
#' @importFrom stats runif
#' @importFrom stats rnorm
#' @export
#' @examples
#' gmtfile <- system.file("extdata", "reactome_subset.gmt", package = "slalom")
#' genesets <- GSEABase::getGmt(gmtfile)
#' data("mesc")
#' model <- newSlalomModel(mesc, genesets, n_hidden = 5, min_genes = 10)
#' model <- init(model)
init <- function(
    object, alpha_priors = NULL, epsilon_priors = NULL,  noise_model = "gauss",
    seed = NULL, pi_prior = NULL, n_hidden = NULL, design = NULL,
    dropFactors = TRUE) {
    ## if Pi priors are supplied, need to redefine object
    if (!is.null(pi_prior)) {
        if (is.null(n_hidden))
            stop("If pi_prior is supplied, n_hidden must be provided.
                 First n_hidden columns of pi_prior are understood as hidden factors. ")
        object <- .newSlalom(object$Y, pi_prior, n_hidden, design)
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
    tmp <- object$Pi_E1
    tmp[tmp < 0.2] <- 0.01
    tmp[tmp > 0.99] <- 0.99
    object$Z_init <- tmp

    ## bits and bobs
    object$onF <- object$nScale # object@nScale
    object$dropFactors <- dropFactors
    # object@isExpressed <- (object@Z[["E1"]] > 0) * 1.0
    # object@numExpressed <- colSums(object@Z[["E1"]] > 0)

    ## initialize parameters using random PCA method
    object <- .initializePCARand_Rcpp(object, seed = seed, saveInit = TRUE)
    ## check validity of object and return
    # validObject(object)
    object
}


.initializePCARand_Rcpp <- function(object, seed = 222, saveInit = FALSE) {
    if (!is.null(seed))
        set.seed(seed)
    # Zstd <- object$Z[["E1"]]
    ## initialize some objects
    Ystd <- object$Y
    unif_vars <- stats::runif(nrow(object$Pi_E1) * ncol(object$Pi_E1))
    Ion <- matrix(unif_vars, nrow = nrow(object$Pi_E1)) < object$Z_init
    object$X_E1 <- matrix(nrow = object$N, ncol = object$K)
    object$X_diagSigmaS <- matrix(1.0 / 2, nrow = object$N, ncol = object$K)
    object$W_E1 <- matrix(nrow = object$G, ncol = object$K)
    ## initialize gamma component of weights
    object$W_gamma0 <- object$Z_init
    object$W_gamma0[object$W_gamma0 <= 0.1] <- 0.1
    object$W_gamma0[object$W_gamma0 >= 0.9] <- 0.9
    object$W_gamma1 <- 1.0 - object$W_gamma0
    ## initialised annotated gene set factors
    for (k in seq_len(object$nAnnotated)) {
        k <- k + object$nKnown
        object$W_E1[, k] <- sqrt(1.0 / object$K) * stats::rnorm(object$G)
        object$X_diagSigmaS[, k] <- 1.0 / 2
        if (sum(Ion[,k]) > 5) {
            ## randomized PCA for speed
            pca <- rsvd::rpca(Ystd[, Ion[, k]], k = 1, retx = TRUE)
            object$X_E1[, k] <- pca$x[, 1] / scale(pca$x[, 1])
        } else {
            object$X_E1[, k] <- stats::runif(object$N)
        }
    }
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
    ## save initial X values if desired
    if (saveInit) {
        object$X_init <- object$X_init
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
#' @param tolerance numeric(1) tolerance to allow between iterations
#' (default 1e-08)
#' @param forceIterations logical(1) should the model be forced to update
#' \code{nIteration} times?
#' @param minIterations integer(1) minimum number of iterations to perform.
#' @param ... generic arguments passed to \code{Rcpp_SlalomModel} method
#'
#' @details Train the model using variational Bayes methods to infer parameters.
#'
#' @docType methods
#' @name train
#' @rdname train
#' @aliases train train,Rcpp_SlalomModel-method
#'
#' @author Davis McCarthy
#' @export
#' @examples
#' gmtfile <- system.file("extdata", "reactome_subset.gmt", package = "slalom")
#' genesets <- GSEABase::getGmt(gmtfile)
#' data("mesc")
#' model <- newSlalomModel(mesc, genesets, n_hidden = 5, min_genes = 10)
#' model <- init(model)
#' model <- train(model, nIterations = 10)
setGeneric("train", function(object, ...) standardGeneric("train"))

#' @name train
#' @rdname train
#' @docType methods
#' @export
setMethod("train", "Rcpp_SlalomModel", function(
    object, nIterations = 5000, tolerance = 1e-08, forceIterations = FALSE,
    minIterations = 700) {
    ## define training parameters in object
    object$tolerance <- tolerance
    object$nIterations <- nIterations
    object$minIterations <- minIterations
    object$forceIterations <- forceIterations
    object$shuffle <- TRUE ## shuffle order in which factors are updated each
    ## for each training iteration

    ## train model
    object$train()
    ## return object
    object
})


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
#'
#' @author Davis McCarthy
#' @export
#' @examples
#' gmtfile <- system.file("extdata", "reactome_subset.gmt", package = "slalom")
#' genesets <- GSEABase::getGmt(gmtfile)
#' data("mesc")
#' model <- newSlalomModel(mesc, genesets, n_hidden = 5, min_genes = 10)
#' model <- init(model)
updateSlalom <- function(object) {
    if (!is(object, "Rcpp_SlalomModel"))
        stop("object must be of class Rcpp_SlalomModel")
    ## train model
    object$update()
    ## return object
    object
}


################################################################################
### Define validity check for SlalomModel class object
#
# setValidity("Rcpp_SlalomModel", function(object) {
#     msg <- NULL
#     valid <- TRUE
#
#     if ( is.null(object@Y) ) {
#         valid <- FALSE
#         msg <- c(msg,
#                  "object must contain an expression matrix at object@Y")
#     }
#     if ( any(is.na(object@Y)) ) {
#         valid <- FALSE
#         msg <- c(msg,
#                  "Expression matrix Y cannot contain NA values")
#     }
#     if ( any(is.na(object@Pi_E1)) ) {
#         valid <- FALSE
#         msg <- c(msg,
#                  "Matrix Pi of observed gene set annotations cannot contain NA values")
#     }
#     if ( (ncol(object@Y) != object@G)) {
#         valid <- FALSE
#         msg <- c(msg,
#                  "Number of genes in expression matrix Y must match object@G")
#     }
#     if ( (nrow(object@Y) != object@N)) {
#         valid <- FALSE
#         msg <- c(msg,
#                  "Number of cells in expression matrix Y must match object@N")
#     }
#     if ( (ncol(object@Pi_E1) != object@K)) {
#         valid <- FALSE
#         msg <- c(msg,
#                  "Number of factors in matrix pi must match object@K")
#     }
#     if ( (nrow(object@Pi_E1) != object@G)) {
#         valid <- FALSE
#         msg <- c(msg,
#                  "Number of factors in matrix pi must match object@G")
#     }
#     if ( (length(object@termNames) != object@K)) {
#         valid <- FALSE
#         msg <- c(msg,
#                  "Number of factors termNames must match object@K")
#     }
#     if ( (length(object@geneNames) != object@G)) {
#         valid <- FALSE
#         msg <- c(msg,
#                  "Number of geneNames must match object@G")
#     }
#     if ( (length(object@cellNames) != object@N)) {
#         valid <- FALSE
#         msg <- c(msg,
#                  "Number of cellNames must match object@N")
#     }
#
#
#     if (valid) TRUE else msg
# })
#
