#' LipidTrend: Analysis and Visualization of Lipid Feature Tendencies
#'
#' @description
#' The LipidTrend package provides tools for analyzing and visualizing trends
#' in lipidomics data. It implements statistical methods for identifying
#' significant changes in lipid abundances across different features,
#' with support for both one-dimensional and two-dimensional analyses.
#'
#' @section Main functions:
#' - \code{\link{analyzeLipidRegion}}: Analyzing lipid trends using permutation
#' tests and Gaussian kernel smoothing
#' - \code{\link{plotRegion1D}}: Visualize one-dimensional lipid trends
#' - \code{\link{plotRegion2D}}: Create two-dimensional visualizations of
#' lipid trends
#'
#' @section Vignettes:
#' See the package vignettes for detailed workflows:
#' \code{vignette('LipidTrend')}
#'
#' @section Installation:
#' To install from Bioconductor, use:
#' \preformatted{
#' if (!requireNamespace('BiocManager', quietly=TRUE))
#'     install.packages('BiocManager')
#' BiocManager::install('LipidTrend')
#' }
#'
#' @name LipidTrend
#' @keywords internal
"_PACKAGE"
