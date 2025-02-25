#' @title Get Result from LipidTrendSE
#' @param object A LipidTrendSE object
#' @return A data frame containing analysis results or NULL
#' @examples
#' data("lipid_se_CL")
#' res_se <- analyzeLipidRegion(
#'     lipid_se=lipid_se_CL, ref_group="sgCtrl", split_chain=FALSE,
#'     chain_col=NULL, radius=3, own_contri=0.5, permute_time=100)
#' # Get complete result
#' results <- getResult(res_se)
#' @export
setGeneric("getResult", function(object) standardGeneric("getResult"))

#' @rdname getResult
#' @exportMethod getResult
setMethod("getResult", "LipidTrendSE", function(object) {
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
#' results <- getEvenChainResult(res_se)
#' @export
setGeneric(
    "getEvenChainResult",
    function(object) standardGeneric("getEvenChainResult"))

#' @rdname getEvenChainResult
#' @aliases getEvenChainResult,LipidTrendSE-method
#' @exportMethod getEvenChainResult
setMethod("getEvenChainResult", "LipidTrendSE", function(object) {
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
#' results <- getOddChainResult(res_se)
#' @export
setGeneric(
    "getOddChainResult",
    function(object) standardGeneric("getOddChainResult"))

#' @rdname getOddChainResult
#' @aliases getOddChainResult,LipidTrendSE-method
#' @exportMethod getOddChainResult
setMethod("getOddChainResult", "LipidTrendSE", function(object) {
    object@odd_chain_result
})


#' @title Get Split Chain Status from LipidTrendSE
#' @param object A LipidTrendSE object
#' @return Logical indicating if analysis was split by chain
#' @keywords internal
setGeneric(
    ".getSplitChainStatus",
    function(object) standardGeneric(".getSplitChainStatus"))

#' @aliases .getSplitChainStatus,LipidTrendSE-method
#' @keywords internal
setMethod(".getSplitChainStatus", "LipidTrendSE", function(object) {
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
        if(.getSplitChainStatus(object)) "Yes" else "No", "\n")

    if (.getSplitChainStatus(object)) {
        even_results <- getEvenChainResult(object)
        if (!is.null(even_results)) {
            cat(
                "Even chain result: ", nrow(even_results),
                " features\n", sep="")
        }
        odd_results <- getOddChainResult(object)
        if (!is.null(odd_results)) {
            cat(
                "Odd chain result: ",
                nrow(odd_results), " features\n", sep="")
        }
    } else {
        result <- getResult(object)
        if (!is.null(result)) {
            cat("Result: ", nrow(result), " features\n", sep="")
        }
    }
})
