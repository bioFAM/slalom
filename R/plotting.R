#### plotting functions for slalom models


#' Plot results of a Slalom model
#'
#' @param object an object of class \code{Rcpp_SlalomModel}
#' @param n_active number of terms (factors) to be plotted (default is 20)
#' @param mad_filter numeric(1), filter factors by this mean absolute deviation
#' to exclude outliers. For large datasets this can be set to 0
#' @param annotated logical(1), should annotated factors be plotted? Default is
#' \code{TRUE}
#' @param unannotated_dense logical(1), should dense unannotated factors be
#' plotted? Default is \code{FALSE}
#' @param unannotated_sparse logical(1), should sparse unannotated factors be
#' plotted? Default is \code{FALSE}
#'
#' @return data.frame with factors ordered by relevance, showing \code{term}
#' (term names), \code{relevance}, \code{type} (factor type: known, annotated
#' or unannotated), \code{n_prior} (number of genes annotated to the gene
#' set/factor), \code{n_gain} (number of genes added/switched on for the
#' factor), \code{n_loss} (number of genes turned off for the factor).
#'
#' @export
#' @examples
#' gmtfile <- system.file("extdata", "reactome_subset.gmt", package = "slalom")
#' genesets <- GSEABase::getGmt(gmtfile)
#' data("mesc")
#' model <- newSlalomModel(mesc, genesets, n_hidden = 5, min_genes = 10)
#' model <- initSlalom(model)
#' model <- trainSlalom(model, nIterations = 10)
#' topTerms(model)
topTerms <- function(
    object, n_active = 20, mad_filter = 0.4, annotated = TRUE,
    unannotated_dense = FALSE, unannotated_sparse = FALSE) {
    if (!methods::is(object, "Rcpp_SlalomModel"))
        stop("object must be of class Rcpp_SlalomModel")

    i_use <- rep(FALSE, object$K)
    if (unannotated_sparse)
        i_use[object$iUnannotatedSparse] <- TRUE
    if (unannotated_dense)
        i_use[object$iUnannotatedDense] <- TRUE
    if (annotated)
        i_use[seq(from = object$nKnown + object$nHidden + 1, to = object$K,
                  by = 1)] <- TRUE

    i_prior <- (object$Pi_E1[, i_use] > 0.5)
    i_posterior <- (object$W_gamma0[, i_use] > 0.5)
    relevance <- (1.0 / object$alpha_E1[i_use])
    terms <- object$termNames[i_use]
    factor_type <- c(rep("known", object$nKnown),
                     rep("unannotated", object$nHidden),
                     rep("annotated", object$K - object$nHidden - object$nKnown))
    factor_type <- factor_type[i_use]
    MAD <- apply(object$X_E1[, i_use], 2, stats::mad)
    R <- (MAD > mad_filter) * relevance

    n_active <- min(sum(R > 0), n_active)

    i_active <- order(R, decreasing = TRUE)[1:n_active]

    df <- data.frame(
        term = terms[i_active],
        relevance = R[i_active],
        type = factor_type[i_active],
        n_prior = colSums(i_prior)[i_active],
        n_gain = colSums(i_posterior & !i_prior)[i_active],
        n_loss = colSums(!i_posterior & i_prior)[i_active]
    )
    df
}


#' Plot results of a Slalom model
#'
#' @param object an object of class \code{Rcpp_SlalomModel}
#' @param n_active number of terms (factors) to be plotted (default is 20)
#' @param mad_filter numeric(1), filter factors by this mean absolute deviation
#' to exclude outliers. For large datasets this can be set to 0
#' @param annotated logical(1), should annotated factors be plotted? Default is
#' \code{TRUE}
#' @param unannotated_dense logical(1), should dense unannotated factors be
#' plotted? Default is \code{FALSE}
#' @param unannotated_sparse logical(1), should sparse unannotated factors be
#' plotted? Default is \code{FALSE}
#'
#' @return a ggplot object
#'
#' @import ggplot2
#' @import grid
#' @export
#' @examples
#' gmtfile <- system.file("extdata", "reactome_subset.gmt", package = "slalom")
#' genesets <- GSEABase::getGmt(gmtfile)
#' data("mesc")
#' model <- newSlalomModel(mesc, genesets, n_hidden = 5, min_genes = 10)
#' model <- initSlalom(model)
#' model <- trainSlalom(model, nIterations = 10)
#' plotRelevance(model)
plotRelevance <- function(
    object, n_active = 20, mad_filter = 0.4, annotated = TRUE,
    unannotated_dense = FALSE, unannotated_sparse = FALSE) {
    if (!methods::is(object, "Rcpp_SlalomModel"))
        stop("object must be of class Rcpp_SlalomModel")

    df <- topTerms(object, n_active, mad_filter, annotated, unannotated_dense,
                   unannotated_sparse)
    df[["term"]] <- factor(df[["term"]], levels = rev(df[["term"]]))
    df[["n_loss"]] <- -df[["n_loss"]]
    p1 <- ggplot(df, aes_string(x = "relevance", y = "term", size = "n_prior")) +
        geom_segment(aes_string(x = 0, xend = "relevance",
                                y = "term", yend = "term"),
                     colour = "gray60", size = 1, show.legend = FALSE) +
        scale_size(range = c(1, 6), name = "Gene set size") +
        xlab("Relevance") +
        ylab("Active pathways") +
        geom_point() +
        theme_classic() + theme(legend.justification = c(1,0),
                                legend.position = c(0.99,0.01))
    p2 <- ggplot(df) +
        geom_segment(aes(x = 0, xend = n_gain,
                                y = term, yend = term),
                     colour = "firebrick", size = 1) +
        geom_point(aes(x = n_gain, y = term, colour = "firebrick"),
                   size = 3) +
        geom_segment(aes(x = 0, xend = n_loss,
                                y = term, yend = term),
                     colour = "dodgerblue", size = 1) +
        geom_point(aes(x = n_loss, y = term, colour = "dodgerblue"),
                   size = 3) +
        scale_colour_manual(name = "Gene set change", guide = "legend",
                            values = c('firebrick' = 'firebrick',
                                       'dodgerblue' = 'dodgerblue'),
                            labels = c('gain','loss')) +
        guides(size = FALSE) +
        xlab("Gene set augmentation") +
        theme_classic() + theme(legend.justification = c(1,0),
                                legend.position = c(0.99,0.01),
                                axis.text.y = element_blank(),
                                axis.title.y = element_blank(),
                                axis.ticks.y = element_blank())

    grid::grid.newpage()
    grid::grid.draw(cbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))
}


#' Plot highest loadings of a factor
#'
#' @param object an object of class \code{Rcpp_SlalomModel}
#' @param term integer(1) or character(1), providing either index for desired
#' term (if an integer) or the term name (if character)
#' @param n_genes integer(1), number of loadings (genes) to show
#'
#' @return a ggplot plot object
#'
#' @export
#' @examples
#' gmtfile <- system.file("extdata", "reactome_subset.gmt", package = "slalom")
#' genesets <- GSEABase::getGmt(gmtfile)
#' data("mesc")
#' model <- newSlalomModel(mesc, genesets, n_hidden = 5, min_genes = 10)
#' model <- initSlalom(model)
#' model <- trainSlalom(model, nIterations = 10)
#' plotLoadings(model, term = 2)
plotLoadings <- function(object, term, n_genes = 10) {
    if (!methods::is(object, "Rcpp_SlalomModel"))
        stop("object must be of class Rcpp_SlalomModel")
    if (is.character(term)) {
        term_match <- object$termNames == term
        if (all(!term_match))
            stop("Provided term not found in object.")
        term <- which(term_match)
    }
    Zchanged <- .getZchanged(object, term)
    W <- object$W_E1[, term]
    Z <- .getZ(object, term)
    gene_labels <- object$geneNames

    Wabs <- abs(W) * abs(Z)
    gene_index <- order(Wabs, decreasing = TRUE)[1:n_genes]

    anno <- rep("in geneset", length(Z))
    anno[Zchanged[gene_index] == 1] <- "gained"

    df <- data.frame(
        gene = gene_labels[gene_index],
        absolute_weight = Wabs[gene_index],
        Annotation = anno[gene_index]
    )
    df[["gene"]] <- factor(df[["gene"]],
                           levels = rev(df[["gene"]]))

    ggplot(df, aes_string(x = "absolute_weight", y = "gene",
                          colour = "Annotation")) +
        geom_segment(aes_string(x = 0, xend = "absolute_weight",
                                y = "gene", yend = "gene"),
                     colour = "gray60", size = 0.5, show.legend = FALSE) +
        scale_size(range = c(1, 6), name = "Gene set size") +
        xlab("Abs. weight") +
        ylab("Genes") +
        geom_point(size = 3) +
        scale_color_manual(values = c("dodgerblue", "firebrick")) +
        theme_classic() + theme(legend.justification = c(1, 0),
                                legend.position = c(0.99,0.01))
}


#' Plot relevance for all terms
#'
#' @param object an object of class \code{Rcpp_SlalomModel}
#' @param terms integer or character vector, providing either indices for
#' desired terms (if an integer) or the term names (if character); default is
#' \code{NULL}, in which case all terms are plotted.
#' @param order_terms logical(1), should factors be ordered by relevance (
#' \code{TRUE}, default), or in the order the come
#' @param mad_filter numeric(1), filter factors by this mean absolute deviation
#' to exclude outliers. For large datasets this can be set close to 0; default
#' is \code{0.2}.
#' @param annotated logical(1), should annotated factors be plotted? Default is
#' \code{TRUE}
#' @param unannotated_dense logical(1), should dense unannotated factors be
#' plotted? Default is \code{TRUE}
#' @param unannotated_sparse logical(1), should sparse unannotated factors be
#' plotted? Default is \code{TRUE}
#'
#' @return a ggplot plot object
#'
#' @export
#' @examples
#' gmtfile <- system.file("extdata", "reactome_subset.gmt", package = "slalom")
#' genesets <- GSEABase::getGmt(gmtfile)
#' data("mesc")
#' model <- newSlalomModel(mesc, genesets, n_hidden = 5, min_genes = 10)
#' model <- initSlalom(model)
#' model <- trainSlalom(model, nIterations = 10)
#' plotTerms(model)
plotTerms <- function(
    object, terms = NULL, order_terms = TRUE, mad_filter = 0.2,
    annotated = TRUE, unannotated_dense = TRUE, unannotated_sparse = FALSE) {
    if (!methods::is(object, "Rcpp_SlalomModel"))
        stop("object must be of class Rcpp_SlalomModel")
    if (is.character(terms)) {
        term_match <- object$termNames == terms
        if (all(!term_match))
            stop("None of the provided terms found in object.")
        terms <- which(term_match)
    }

    df <- topTerms(object, n_active = object$K, mad_filter, annotated,
                   unannotated_dense, unannotated_sparse)
    if (order_terms)
        df[["term"]] <- factor(df[["term"]], levels = rev(df[["term"]]))
    else
        df[["term"]] <- factor(df[["term"]], levels = rev(sort(df[["term"]])))

    ggplot(df, aes_string(x = "relevance", y = "term",
                          colour = "type")) +
        geom_segment(aes_string(x = 0, xend = "relevance",
                                y = "term", yend = "term"),
                     colour = "gray80", size = 0.5, show.legend = FALSE) +
        scale_size(range = c(1, 4), name = "Gene set size") +
        xlab("Relevance") +
        ylab("Terms") +
        geom_point(size = 3) +
        scale_color_manual(values = c("dodgerblue", "firebrick"),
                           name = "Term type") +
        theme_classic() + theme(legend.justification = c(1, 0),
                                legend.position = c(0.99,0.01))
}



## Get posterior of Z (Bernourlli part part of spike-and-slab prior) :math:`Q(Z)`
## Args: term        (str): optional list of terms for which weights are returned. Default None=all terms.
.getZ <- function(object, terms = NULL) {
    if (is.null(terms))
        object$W_gamma0
    else
        object$W_gamma0[, terms]
}


## get matrix indicating whether the posterior distribution has changed for individual terms/genes
## Args: terms        (str): optional list of terms for which weights are returned. Default None=all terms.
## Rv: matrix [0,-1,1]: 0: no change, -1: loss, +1: gain
.getZchanged <- function(object, terms = NULL, threshold = 0.5) {
    if (is.null(terms)) terms <- seq_len(object$K)
    Z <- .getZ(object, terms)
    Pi <- object$Pi_E1[, terms]
    I <- matrix(0, nrow = object$G, ncol = length(terms))
    Igain <- (Z > threshold) & (Pi < threshold)
    Iloss <- (Z < threshold) & (Pi > threshold)
    I[Igain] <- 1
    I[Iloss] <- -1
    I
}





