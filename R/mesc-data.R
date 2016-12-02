#' @name mesc
#' @title A single-cell expression dataset to demonstrate
#' capabilities of slalom from mouse embryonic stem cells (mESCs)
#' @description This data set consists of an \code{\link[scater]{SCESet}} object
#' with log2-counts-per-million expression values for 3635 genes for 182 cells.
#' They are from a real experiment, studying cell cycle in mouse embryonic stem
#' cells (mESCs). See Buettner et al (Nat. Biotech., 2015) for details. d.
#' @return NULL, but makes aavailable an SCESet object containing expression data
#' @docType data
#' @usage mesc
#' @format an SCESet instance, 1 row per gene.
#' @source EMBL-EBI, Hinxton, UK
#' @references
#' Buettner F, Natarajan KN, Paolo Casale F, Proserpio V, Scialdone A, Theis FJ, et al. Computational analysis of cell-to-cell heterogeneity in single-cell RNA-sequencing data reveals hidden subpopulations of cells. Nat Biotechnol. Nature Publishing Group; 2015;33: 155–160.
#' @author Davis McCarthy, Florian Buettner, 2016-12-02
NULL