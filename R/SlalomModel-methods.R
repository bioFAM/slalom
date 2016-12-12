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
#'
#' exprsfile <- system.file("extdata", "mesc.csv", package = "slalom")
#' mesc_mat <- as.matrix(read.csv(exprsfile))
#' sce <- scater::newSCESet(exprsData = mesc_mat, logExprsOffset = 1,
#'                  lowerDetectionLimit = 0)
#' # model2 <- newSlalomModel(mesc_mat, genesets, n_hidden = 5, min_genes = 10)
#'
newSlalomModel <- function(object, genesets, n_hidden = 5, prune_genes = TRUE,
                           min_genes = 15, design = NULL) {
    n_known <- ifelse(is.null(design), 0, ncol(design))
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
    ## set some attributes that we need frequently for the updates, inculuding
    ## number and idx of hidden (unannotated) terms
    out <- new("SlalomModel",
        K = n_known + ncol(I),
        N = nrow(Y),
        G = nrow(I),
        nAnnotated = ncol(I) - n_hidden,
        nHidden = n_hidden,
        nKnown = n_known,
        nScale = 100,
        Y = Y[, retained_genes],
        pi = I * 1L,
        termNames = colnames(I),
        geneNames = retained_genes,
        cellNames = scater::cellNames(object),
        iUnannotatedDense = seq(from = (ncol(I) - n_hidden + 1), to = ncol(I))
    )
    if (!is.null(design)) {
        if (nrow(design) != out@N)
            stop("Number of rows of design does not match number of cells.")
        else
            out@Known <- design
    }
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
    if ( any(is.na(object@pi)) ) {
        valid <- FALSE
        msg <- c(msg,
                 "Matrix pi of observed gene set annotations cannot contain NA values")
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
    if ( (ncol(object@pi) != object@K)) {
        valid <- FALSE
        msg <- c(msg,
                 "Number of factors in matrix pi must match object@K")
    }
    if ( (nrow(object@pi) != object@G)) {
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


################################################################################
### Initialize a SlalomModel class object

#' Initialize a SlalomModel object
#'
#' Initialize a SlalomModel with sensible starting values for parameters before
#' training the model.
#'
#' @param object a \code{SlalomModel} object
#' @param alpha_priors numeric(2) giving alpha and beta hyperparameters for a
#' gamma prior distribution for alpha parameters (precision of factor weights)
#' @param epsilon_priors numeric(2) giving alpha and beta hyperparameters for a
#' gamma prior distribution for noise precision parameters
#' @param noise character(1) defining noise model, defaults to "gauss" for
#' Gaussian noise model
#' @param ... generic arguments passed to
#'
#' @details It is strongly recommended to use \code{\link{newSlalomModel}} to
#' create the \code{\link{SlalomModel}} object prior to applying
#' \code{initialize}.
#'
#' @docType methods
#' @name init
#' @rdname init
#' @aliases init init,SlalomModel-method
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
setGeneric("init", function(object, ...) standardGeneric("init"))

#' @name init
#' @rdname init
#' @docType methods
#' @export
setMethod("init", "SlalomModel", function(
    object, alpha_priors = NULL, epsilon_priors = NULL,  noise = "gauss") {
    ## define priors for alpha and epsilon
    if (is.null(object@alpha[["priors"]]))
        object@alpha[["priors"]] <- c(1e-03, 1e-03)
    if (is.null(object@epsilon[["priors"]]))
        object@epsilon[["priors"]] <- c(1e-03, 1e-03)

    ## pi is likelihood of link for genes x factors, defined in object from
    ## annotated gene sets with newSlalomModel
    ## define the number of genes that are "on" for each factor
    object@nOn <- colSums(object@pi > 0.5)

    ## define initial expected values for indicator variables Z with observed
    ## annotations
    object@Z[["E1"]] <- object@pi
    object@Z[["initZ"]] <- object@pi
    object@Z[["initZ"]][object@Z[["initZ"]] < 0.2] <- 0.01

    ## bits and bobs
    # object@onF <- object@nScale # object@nScale
    # object@isExpressed <- (object@Z[["E1"]] > 0) * 1.0
    # object@numExpressed <- colSums(object@Z[["E1"]] > 0)

    ## initialize parameters using random PCA method
    object <- .initializePCARand(object)
    ## check validity of object and return
    validObject(object)
    object
})


.initializePCARand <- function(object, seed = 222, saveInit = FALSE) {
    set.seed(seed)
    # Zstd <- object@Z[["E1"]]
    ## initialize some objects
    Ystd <- object@Y
    unif_vars <- stats::runif(nrow(object@pi) * ncol(object@pi))
    Ion <- matrix(unif_vars, nrow = nrow(object@pi)) < object@Z[["initZ"]]
    object@X[["E1"]] <- matrix(nrow = object@N, ncol = object@K)
    object@X[["diagSigmaS"]] <- matrix(1.0 / 2, nrow = object@N, ncol = object@K)
    object@W[["E1"]] <- matrix(nrow = object@G, ncol = object@K)
    ## initialize gamma component of weights
    object@W[["gamma0"]] <- object@Z[["initZ"]]
    object@W[["gamma0"]][object@W[["gamma0"]] <= 0.1] <- 0.1
    object@W[["gamma0"]][object@W[["gamma0"]] >= 0.9] <- 0.9
    object@W[["gamma1"]] <- 1.0 - object@W[["gamma0"]]
    ## initialised annotated gene set factors
    for (k in seq_len(object@nAnnotated)) {
        k <- k + object@nKnown
        if (sum(Ion[,k]) > 5) {
            ## randomized PCA for speed
            pca <- rsvd::rpca(Ystd[, Ion[, k]], k = 1, retx = TRUE)
            object@X[["E1"]][, k] <- pca$x[, 1] / scale(pca$x[, 1])
        } else {
            object@X[["E1"]][, k] <- stats::runif(object@N)
            object@W[["E1"]][, k] <- sqrt(1.0 / object@K) * stats::rnorm(object@G)
            object@X[["diagSigmaS"]][, k] <- 1.0 / 2
        }
    }
    ## initialise known factors
    if (object@nKnown > 0L) {
        for (k in seq_len(object@nKnown)) {
            object@W[["E1"]][, k] <- sqrt(1.0 / object@K) * stats::rnorm(object@G)
            object@X[["diagSigmaS"]][, k] <- 1.0 / 2
        }
        object@X[["E1"]][, seq_len(object@nKnown)] <- object@Known
    }
    ## initialise hidden (latent) factors
    if (object@nHidden > 0L) {
        for (iL in object@iUnannotatedDense)
            object@X[["E1"]][,iL] <- stats::rnorm(object@N)
    }
    ## save initial X values - currently ignored
    if (saveInit) {
        object@X[["initX"]] <- object@X[["E1"]]
    }
    object
}


# .initializePCA <- function(object, seed = 222) {
#     set.seed(seed)
#     #pca initialisation
#     Ion <- matrix(stats::runif(nrow(object@pi) * ncol(object@pi))) < object@pi
#     object@W[["gamma"]] <- object@pi
#     for (k in seq_len(object@K)) {
#         sv <- svd(object@Z[["E1"]][, Ion[:, k]])
#         s0 <- sv$u[, 1]
#         w0 <- t(sv$d %*% sv$v)[, 1]
#         v <- scale(s0)
#         s0 <- s0 / v
#         w0 <- w0 * v
#         object@X[["E1"]][, k] <- s0
#         object@W[["E1"]][Ion[, k], k] <- w0
#         object@W[["E1"]][!Ion[, k], k] <- (object@sigmaOff *
#                                                object@W[["E1"]][!Ion[, k], k])
#         object@X[["diagSigmaS"]][, k] <- 1.0 / 2
#     }
#     object
# }

#
# .initializeGreedy <- function() {
#     object@X[["E1"]] = random.randn(self._N,object@K)
#     self.W.E1 = random.randn(self._D,object@K)
#     Ion = (self.Pi>0.5)
#     self.W.E1[~Ion]*= self.sigmaOff
#     for k in range(Ion.shape[1]):
#         self.W.E1[Ion[:,k]]*=self.sigmaOn[k]
# }
#
#
# .initializePrior <- function() {
#     Ion = random.rand(self.Pi.shape[0],self.Pi.shape[1])<self.Pi
#     self.W.E1[~Ion]*=self.sigmaOff
#     for k in range(Ion.shape[1]):
#         self.W.E1[Ion[:,k],k]*=self.sigmaOn[k]
# }
#
#
# .initializeOn <- function() {
#     for k in range(Ion.shape[1]):
#         self.W.E1[:,k]*=self.sigmaOn[k]
# }
#
#
# .initializeRandom <- function() {
#     for k in range(self.Pi.shape[1]):
#         self.S.diagSigmaS[:,k] = 1./2
#         object@X[["E1"]][:,k] = SP.randn(self._N)
#         self.W.E1 = SP.randn(self._D, self.Pi.shape[1])
#         object@W[["gamma"]][:,:,0] = self.Pi
#         object@W[["gamma"]][:,:,0][object@W[["gamma"]][:,:,0]<=.2] = .1
#         object@W[["gamma"]][:,:,0][object@W[["gamma"]][:,:,0]>=.8] = .9
#         if self.nKnown>0:
#             for k in SP.arange(self.nKnown):
#             self.W.E1[:,k] = SP.sqrt(1./object@K)*SP.randn(self._D)
#         self.S.diagSigmaS[:,k] = 1./2
#         object@X[["E1"]][:,SP.arange(self.nKnown)] =  self.Known
#         if self.saveInit==True:
#             self.initS = object@X[["E1"]].copy()
#
# }
#
# .initializeData <- function() {
#     assert ('S' in list(init_factors.keys()))
#     assert ('W' in list(init_factors.keys()))
#     #            Ion = init_factors['Ion']
#     Sinit = init_factors['S']
#     Winit = init_factors['W']
#     object@W[["gamma"]][:,:,0] = self.Pi
#     object@W[["gamma"]][:,:,0][object@W[["gamma"]][:,:,0]<=.2] = .1
#     object@W[["gamma"]][:,:,0][object@W[["gamma"]][:,:,0]>=.8] = .9
#     for k in range(object@K):
#         object@X[["E1"]][:,k] = Sinit[:,k]
#     self.W.E1[:,k] = Winit[:,k]
#     self.S.diagSigmaS[:,k] = 1./2
# }
#


