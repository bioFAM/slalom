# Define the SlalomModel class for the package

#' The "Slalom Model" (SlalomModel)  class
#'
#' S4 class and the main class used by scater to hold single cell expression
#' data. SCESet extends the basic Bioconductor ExpressionSet class.
#'
#' This class is initialized from a matrix of expression values.
#'
#' Methods that operate on SCESet objects constitute the basic scater workflow.
#'
#'@section Slots:
#'  \describe{
#'    \item{\code{K}:}{Scalar of class \code{"numeric"}, indicating total number
#'     of factors.}
#'    \item{\code{N}:}{Scalar of class \code{"numeric"}, indicating number of
#'    cells.}
#'    \item{\code{G}:}{Scalar of class \code{"numeric"}, indicating number of
#'    genes.}
#'    \item{\code{nScale}:}{Scalar of class \code{"numeric"}, relative number of
#'     cells to which the observed annotation data (gene sets) are scaled.}
#'    \item{\code{nAnnotated}:}{Scalar of class \code{"numeric"}, number of
#'    factors relating to annotated gene sets to model.}
#'    \item{\code{nHidden}:}{Scalar of class \code{"numeric"}, number of hidden
#'    (latent) factors to model.}
#'    \item{\code{nKnown}:}{Scalar of class \code{"numeric"}, number of known
#'    factors (covariates).}
#'    \item{\code{noiseModel}:}{\code{"character"} scalar defining noise model
#'    used by the model (default: "gauss" for Gaussian noise model.}
#'    \item{\code{iUnannotatedDense}:}{\code{"integer"} vector giving indices for
#'    dense unannotated (hidden) factors.}
#'    \item{\code{iUnannotatedSparse}:}{\code{"integer"} vector giving indices for
#'    sparse unannotated (hidden) factors.}
#'    \item{\code{nOn}:}{Vector of class \code{"numeric"}, number of genes that
#'    are "on" for each factor (as annotated).}
#'    \item{\code{termNames}:}{\code{"character"} vector giving names for the gene
#'    sets and factors.}
#'    \item{\code{geneNames}:}{\code{"character"} vector giving gene names.}
#'    \item{\code{cellNames}:}{\code{"character"} vector giving cell names.}
#'    \item{\code{X}:}{\code{"list"} of values and parameters for factor states.}
#'    \item{\code{W}:}{\code{"list"} of values and parameters for factor weights.}
#'    \item{\code{Z}:}{\code{"list"} of values and parameters for activation
#'    variable of whether factor k has a regulatory effect on gene g.}
#'    \item{\code{alpha}:}{\code{"list"} of values and parameters for factor
#'    precisions.}
#'    \item{\code{epsilon}:}{\code{"list"} of values and parameters for residual
#'    precisions.}
#'    \item{\code{pi}:}{\code{"matrix"} of size G x K with each entry being the
#'    prior probability for a gene g being active for factor k.}
#'    \item{\code{I}:}{\code{"matrix"} of size G x K of observed annotation data
#'     with each entry being the indicator that gene g is annotated to factor k.}
#'    \item{\code{Known}:}{design \code{"matrix"} defining covariates to fit in
#'    the model ("known factors").}
#'    \item{\code{Y}:}{\code{"matrix"} of size N x G with each entry being the
#'    observed expression value (normalized, log2-scale) for gene g in cell n.}
#'    \item{\code{pseudo_Y}:}{\code{"matrix"} of size N x G with each entry
#'    being the pseudoexpression value (normalized, log2-scale) for gene g in
#'    cell n.}
#'}
#' @name SlalomModel
#' @rdname SlalomModel
#' @aliases SlalomModel-class
#' @exportClass SlalomModel
setClass("SlalomModel",
         slots = c(K = "numeric",
                   N = "numeric",
                   G = "numeric",
                   nScale = "numeric",
                   nAnnotated = "numeric",
                   nHidden = "numeric",
                   nKnown = "numeric",
                   noiseModel = "character",
                   iUnannotatedDense = "integer",
                   iUnannotatedSparse = "integer",
                   nOn = "numeric",
                   termNames = "character",
                   geneNames = "character",
                   cellNames = "character",
                   X = "list",
                   W = "list",
                   Z = "list",
                   alpha = "list",
                   epsilon = "list",
                   pi = "matrix",
                   I = "matrix",
                   Known = "matrix",
                   Y = "matrix",
                   pseudo_Y = "matrix")
)


