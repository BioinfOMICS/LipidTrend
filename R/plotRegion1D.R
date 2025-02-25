#' @title Plot region trends for one-dimensional lipid feature
#' @description Visualize analysis results for one-dimensional features to
#' identify lipid feature tendencies. The plots highlight regions of significant
#' differences between case and control samples: blue areas indicate
#' significantly lower abundance in case samples, while red areas indicate
#' significantly higher abundance.
#' @param object A LipidTrendSE object containing analysis results.
#' @param p_cutoff Numeric. Significance level. The range is 0 to 1.
#'   Default is \code{0.05}.
#' @return For split chain analyses, a list with two elements:
#'   \item{even_result}{Plot for even chain features, or NULL if none exist}
#'   \item{odd_result}{Plot for odd chain features, or NULL if none exist}
#'   For non-split analyses, a single plot.
#' @examples
#' data("lipid_se_CL")
#' res_se <- analyzeLipidRegion(
#'     lipid_se=lipid_se_CL, ref_group="sgCtrl", split_chain=FALSE,
#'     chain_col=NULL, radius=2, own_contri=0.5, permute_time=100,
#'     abun_weight=TRUE)
#' plot <- plotRegion1D(res_se, p_cutoff=0.05)
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line geom_point labs
#'   scale_color_manual theme_bw
#' @importFrom rlang sym !!
#' @export
#' @seealso
#' \code{\link{analyzeLipidRegion}} for generating the input LipidTrendSE object
setGeneric(
    "plotRegion1D",
    function(object, p_cutoff = 0.05) standardGeneric("plotRegion1D"))

#' @rdname plotRegion1D
#' @aliases plotRegion1D,LipidTrendSE-method
#' @exportMethod plotRegion1D
setMethod("plotRegion1D", "LipidTrendSE", function(object, p_cutoff=0.05) {
    if (!is.numeric(p_cutoff) || p_cutoff < 0 || p_cutoff > 1) {
        stop("p_cutoff must be a numeric value between 0 and 1")
    }
    if (.getSplitChainStatus(object)) {
        even_results <- getEvenChainResult(object)
        odd_results <- getOddChainResult(object)
        plots <- list(
            even_result=if(
                !is.null(even_results)) .trendPlot1D(
                    even_results, p_cutoff) else NULL,
            odd_result=if(
                !is.null(odd_results)) .trendPlot1D(
                    odd_results, p_cutoff) else NULL
        )
        return(plot1D=plots)
    } else {
        results <- getResult(object)
        plot <- if(
            !is.null(results)) .trendPlot1D(results, p_cutoff) else NULL
        return(plot1D=plot)
    }
})

#' Create 1D result plot
#' @description Creates a ggplot visualization showing trends between case and
#'   control groups with highlighted significant regions.
#' @param res A data.frame containing columns for coordinates, abundance values,
#'   direction, and adjusted p-values
#' @param p_cutoff Numeric value for significance cutoff (default=0.05)
#' @return A ggplot2 object
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line geom_point labs
#'   scale_color_manual theme_bw
#' @importFrom rlang sym !!
#' @keywords internal
.trendPlot1D <- function(res, p_cutoff) {
    result <- ggplot2::ggplot(
        res, ggplot2::aes(x = !!rlang::sym(colnames(res)[1])))
    num.neg.region <- sum(
        res$direction == "-" & res$smoothing.pval.BH<p_cutoff)
    if (num.neg.region > 0) {
        result <- result + ggplot2::geom_ribbon(
            ggplot2::aes(ymin=ifelse(
                smoothing.pval.BH<p_cutoff & avg.expr.ctrl>=avg.expr.case,
                avg.expr.case, NA), ymax=avg.expr.ctrl), fill="blue", alpha=0.3)
    }
    num.pos.region <- sum(
        res$direction == "+" & res$smoothing.pval.BH<p_cutoff)
    if (num.pos.region>0) {
        result <- result + ggplot2::geom_ribbon(
            ggplot2::aes(ymin=ifelse(
                smoothing.pval.BH<p_cutoff & avg.expr.case>=avg.expr.ctrl,
                avg.expr.ctrl, NA), ymax=avg.expr.case), fill="red", alpha=0.3)
    }
    result <- result +
        ggplot2::geom_line(ggplot2::aes(y=avg.expr.ctrl, color="Control")) +
        ggplot2::geom_line(ggplot2::aes(y=avg.expr.case, color="Case")) +
        ggplot2::geom_point(ggplot2::aes(y=avg.expr.ctrl, color="Control")) +
        ggplot2::geom_point(ggplot2::aes(y=avg.expr.case, color="Case")) +
        ggplot2::labs(color="Group", y="Lipid Abundance") +
        ggplot2::scale_color_manual(
            breaks=c("Control", "Case"), values=c("blue", "red")) +
        ggplot2::theme_bw()
    return(result)
}
