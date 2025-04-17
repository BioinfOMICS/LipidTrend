#' @title Get Result from LipidTrendSE
#' @param object A LipidTrendSE object
#' @return A data frame containing analysis results or NULL
#' @examples
#' data("lipid_se_CL")
#' res_se <- analyzeLipidRegion(
#'     lipid_se=lipid_se_CL, ref_group="sgCtrl", split_chain=FALSE,
#'     chain_col=NULL, radius=3, own_contri=0.5, permute_time=100)
#' # Get complete result
#' results <- result(res_se)
#' @export
setGeneric("result", function(object) standardGeneric("result"))

#' @rdname result
#' @exportMethod result
setMethod("result", "LipidTrendSE", function(object) {
    object@result
})

#' @title Get Even Chain Result from LipidTrendSE
#' @param object A LipidTrendSE object
#' @return A data frame containing even chain result or NULL
#' @examples
#' data("lipid_se_CL")
#' sub <- lipid_se_CL[seq_len(10), ]
#' res_se <- analyzeLipidRegion(
#'     lipid_se=sub, ref_group="sgCtrl", split_chain=TRUE,
#'     chain_col="chain", radius=3, own_contri=0.5, permute_time=100)
#' # Get complete result summary
#' results <- even_chain_result(res_se)
#' @export
setGeneric(
    "even_chain_result",
    function(object) standardGeneric("even_chain_result"))

#' @rdname even_chain_result
#' @aliases even_chain_result,LipidTrendSE-method
#' @exportMethod even_chain_result
setMethod("even_chain_result", "LipidTrendSE", function(object) {
    object@even_chain_result
})

#' @title Get Odd Chain Result from LipidTrendSE
#' @param object A LipidTrendSE object
#' @return A data frame containing odd chain result or NULL
#' @examples
#' data("lipid_se_CL")
#' res_se <- analyzeLipidRegion(
#'     lipid_se=lipid_se_CL, ref_group="sgCtrl", split_chain=TRUE,
#'     chain_col="chain", radius=3, own_contri=0.5, permute_time=100)
#' # Get complete result summary
#' results <- odd_chain_result(res_se)
#' @export
setGeneric(
    "odd_chain_result",
    function(object) standardGeneric("odd_chain_result"))

#' @rdname odd_chain_result
#' @aliases odd_chain_result,LipidTrendSE-method
#' @exportMethod odd_chain_result
setMethod("odd_chain_result", "LipidTrendSE", function(object) {
    object@odd_chain_result
})


#' @title Get Split Chain Status from LipidTrendSE
#' @param object A LipidTrendSE object
#' @return Logical indicating if analysis was split by chain
#' @keywords internal
setGeneric(
    ".split_chain",
    function(object) standardGeneric(".split_chain"))

#' @aliases .split_chain,LipidTrendSE-method
#' @keywords internal
setMethod(".split_chain", "LipidTrendSE", function(object) {
    object@split_chain
})

#' @title Show method for LipidTrendSE objects
#' @param object A LipidTrendSE object
#' @return LipidTrendSE object information
#' @importMethodsFrom SummarizedExperiment show
#' @aliases show,LipidTrendSE-method
#' @exportMethod show
setMethod("show", "LipidTrendSE", function(object) {
    callNextMethod()
    cat("\nLipidTrend Results:\n")
    cat("------------------------\n")
    cat("Split chain analysis:",
        if(.split_chain(object)) "Yes" else "No", "\n")

    if (.split_chain(object)) {
        even_results <- even_chain_result(object)
        if (!is.null(even_results)) {
            cat(
                "Even chain result: ", nrow(even_results),
                " features\n", sep="")
        }
        odd_results <- odd_chain_result(object)
        if (!is.null(odd_results)) {
            cat(
                "Odd chain result: ",
                nrow(odd_results), " features\n", sep="")
        }
    } else {
        result <- result(object)
        if (!is.null(result)) {
            cat("Result: ", nrow(result), " features\n", sep="")
        }
    }
})
