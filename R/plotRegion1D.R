#' @title Plot region trends for one-dimensional lipid feature
#' @description Visualize analysis results for one-dimensional features to
#' identify lipid feature tendencies. The plots highlight regions of significant
#' differences between case and control samples: blue areas indicate
#' significantly lower abundance in case samples, while red areas indicate
#' significantly higher abundance.
#' @param object A LipidTrendSE object containing analysis results.
#' @param p_cutoff Numeric. Significance level. The range is 0 to 1.
#'   Default is \code{0.05}.
#' @param y_scale Character. Choose one of the y-axis scales: "identity",
#' "log2", "log10", or "sqrt". Default is \code{"identity"}.
#' @return For split chain analyses, a list with two elements:
#'   \item{even_result}{Plot for even chain features, or NULL if none exist}
#'   \item{odd_result}{Plot for odd chain features, or NULL if none exist}
#'   For non-split analyses, a single plot.
#' @examples
#' data("lipid_se_CL")
#' res_se <- analyzeLipidRegion(
#'     lipid_se=lipid_se_CL, ref_group="sgCtrl", split_chain=FALSE,
#'     chain_col=NULL, radius=2, own_contri=0.5, permute_time=100,
#'     abund_weight=TRUE)
#' plot <- plotRegion1D(res_se, p_cutoff=0.05, y_scale='identity')
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line geom_point labs
#'   scale_color_manual theme_bw scale_x_continuous scale_y_continuous
#'   scale_fill_manual element_text element_line element_blank
#' @importFrom rlang sym !!
#' @export
#' @seealso
#' \code{\link{analyzeLipidRegion}} for generating the input LipidTrendSE object
setGeneric(
    "plotRegion1D",
    function(
        object, p_cutoff=0.05, y_scale='identity')
        standardGeneric("plotRegion1D"))

#' @rdname plotRegion1D
#' @aliases plotRegion1D,LipidTrendSE-method
#' @exportMethod plotRegion1D
setMethod("plotRegion1D", "LipidTrendSE", function(
        object, p_cutoff=0.05, y_scale='identity') {
    if (!is.numeric(p_cutoff) || p_cutoff < 0 || p_cutoff > 1) {
        stop("p_cutoff must be a numeric value between 0 and 1")
    }
    if (is.null(y_scale) ||
        !(y_scale %in% c('identity', 'log2', 'log10', 'sqrt')) ) {
        stop("y_scale must be one of identity, log2, log10, or sqrt")
    }
    if (.split_chain(object)) {
        even_results <- even_chain_result(object)
        odd_results <- odd_chain_result(object)
        plots <- list(
            even_result=if (!is.null(even_results)) {
                .trendPlot1D(even_results, p_cutoff, y_scale) +
                    scale_x_continuous(breaks=seq(0,150,2))
            } else NULL,
            odd_result=if (!is.null(odd_results)) {
                .trendPlot1D(odd_results, p_cutoff, y_scale) +
                    scale_x_continuous(breaks=seq(1,150,2))
            } else NULL
        )
        return(plot1D=plots)
    } else {
        results <- result(object)
        plot <- if (!is.null(results)) {
            .trendPlot1D(results, p_cutoff, y_scale) +
                scale_x_continuous(breaks=seq(0,150,1))
        } else NULL
        return(plot1D=plot)
    }
})

#' Create 1D result plot
#' @description Creates a ggplot visualization showing trends between case and
#'   control groups with highlighted significant regions.
#' @param res A data.frame containing columns for coordinates, abundance values,
#'   direction, and adjusted p-values
#' @param p_cutoff Numeric value for significance cutoff.
#' Default is \code{0.05}
#' @param y_scale Character. Choose one of the y-axis scales: "identity",
#' "log2", "log10", or "sqrt". Default is \code{"identity"}
#' @return A ggplot2 object
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line geom_point labs
#'   scale_color_manual theme_bw scale_x_continuous scale_y_continuous
#'   scale_fill_manual element_text element_line element_blank
#' @importFrom rlang sym !!
#' @keywords internal
.trendPlot1D <- function(res, p_cutoff, y_scale) {
    result <- ggplot(res, aes(x = !!rlang::sym(colnames(res)[1])))
    num.neg.region <- sum(
        res$direction == "-" & res$smoothing.pval.BH < p_cutoff)
    if (num.neg.region > 0) {
        result <- result + geom_ribbon(
            aes(ymin=ifelse(
                smoothing.pval.BH < p_cutoff & avg.abund.ctrl >= avg.abund.case,
                avg.abund.case, NA),
                ymax=avg.abund.ctrl, fill="Low"), alpha=0.3)
    }
    num.pos.region <- sum(
        res$direction == "+" & res$smoothing.pval.BH < p_cutoff)
    if (num.pos.region > 0) {
        result <- result + geom_ribbon(
            aes(ymin=ifelse(
                smoothing.pval.BH < p_cutoff & avg.abund.case >= avg.abund.ctrl,
                avg.abund.ctrl, NA), ymax=avg.abund.case, fill="High"),
            alpha=0.3)
    }
    result <- result +
        geom_line(aes(y=avg.abund.ctrl, color="Control"), linewidth=0.8) +
        geom_line(aes(y=avg.abund.case, color="Case"), linewidth=0.8) +
        geom_point(aes(
            y=avg.abund.ctrl, color="Control"), size=2.3, shape=16) +
        geom_point(aes(y=avg.abund.case, color="Case"), size=3, shape=18) +
        scale_color_manual(
            breaks=c("Case", "Control"), values=c("#FF5151", "#4169E1")) +
        scale_fill_manual(
            breaks = c("High", "Low"),
            values=c("High"="#FF5151", "Low"="#4169E1")) +
        labs(color="Group", fill="Trend", y="Lipid Abundance")
    result <- result +
        (switch(y_scale, "log2"=scale_y_continuous(trans = "log2"),
                "log10"=scale_y_continuous(trans = "log10"),
                "sqrt"=scale_y_continuous(trans = "sqrt"),
                "identity"=scale_y_continuous()))
    result <- result + theme_bw() +
        theme(
            axis.title=element_text(size=15), axis.text=element_text(size=14),
            legend.title=element_text(size=12),
            legend.text=element_text(size=12), panel.grid.minor=element_blank(),
            panel.grid.major=element_line(color="grey95"))
    return(result)
}
