#' @title Conduct statistical to analyze lipid features tendencies
#' @description Calculate p-value of lipidomic features by permutation test.
#' The test of region statistics smoothing with Gaussian kernel integrates
#' neighbor's information and provides a stable testing result under small
#' sample size.
#' @param lipid_se A SummarizedExperiment object. Must contain following data:
#' \describe{
#'      \item{Assay}{A matrix containing lipid abundance data, where rows
#'      represent lipids and columns represent samples.}
#'      \item{RowData}{A data frame of lipid features (e.g., double bond count,
#'      chain length), with rows as lipids and columns as lipid features (
#'      limited to 1 or 2 columns). The order of lipids must match the abundance
#'      data. If RowData contains one column, a one-dimensional analysis will
#'      be performed. If RowData includes two columns, a two-dimensional
#'      analysis will be conducted.}
#'      \item{ColData}{A data frame containing group information, where rows
#'      represent sample names and columns must include sample name, label name,
#'      and group, arranged accordingly.}
#' }
#' @param ref_group Character. Group name of the reference group. It must be
#' one of the group names in the colData group column.
#' @param split_chain Logical. If TRUE the results will split to shown by odd
#' and even chain. Default is \code{FALSE}.
#' @param chain_col Character. The column name of chain length. Set to NULL if
#' split_chain is FALSE. Default is \code{'NULL'}.
#' @param radius Numeric. Distance of neighbor to be included.
#' Default is \code{3}.
#' @param own_contri Numeric. Proportion of own contribution.
#' Default is \code{0.5}. To not over-emphasize the neighbor information,
#' we suggest to choose a value from 0.5 to 1.
#' @param test string. Choose statistic test from "t.test" and "Wilcoxon".
#' Default is \code{'t.test'}.
#' @param abun_weight Logical. Consider the average abundance as the weight in
#' the test statistic. Default is \code{TRUE}.
#' @param permute_time Integer. Permutation times. Default is \code{100000}.
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
#' @importFrom SummarizedExperiment assay rowData colData
#' @return A LipidTrendSE object containing lipidomic feature testing result.
#' @examples
#' data("lipid_se_CL")
#' res_se <- analyzeLipidRegion(
#'     lipid_se=lipid_se_CL, ref_group="sgCtrl", split_chain=FALSE,
#'     chain_col=NULL, radius=3, own_contri=0.5, permute_time=100)
#' @export
#' @seealso
#' \code{\link{plotRegion1D}} for one-dimensional visualization
#' \code{\link{plotRegion2D}} for two-dimensional visualization
analyzeLipidRegion <- function(
        lipid_se, ref_group, split_chain=FALSE, chain_col=NULL, radius=3,
        own_contri=0.5, test="t.test", abun_weight=TRUE, permute_time=100000){
    if (isFALSE(inherits(lipid_se, "SummarizedExperiment")) ) {
        stop("Invalid lipid_se input.")
    }
    if (!is.numeric(radius) || radius <= 0) {
        stop("radius must be a positive numeric value")
    }
    if (!is.numeric(own_contri) || own_contri < 0 || own_contri > 1) {
        stop("own_contri must be between 0 and 1")
    }
    X <- assay(lipid_se) %>% t()
    X.info <- rowData(lipid_se) %>% as.data.frame()
    group_info <- colData(lipid_se) %>% as.data.frame()
    if (is.null(ref_group) || !(ref_group %in% group_info$group)) {
        stop("ref_group must be one of the groups in colData group column")
    }
    group_info %<>% dplyr::mutate(
        lipidTrend_group=ifelse(group_info$group==ref_group, 0, 1))
    group <- group_info$lipidTrend_group
    if (isTRUE(split_chain)) {
        if (!is.null(chain_col) && chain_col %in% colnames(X.info)) {
            chain_parity <- X.info[[chain_col]] %% 2
            chain_groups <- list(even=chain_parity == 0, odd=chain_parity == 1)
            output.df <- lapply(names(chain_groups), function(chain_type) {
                sub_feature <- chain_groups[[chain_type]]
                if (sum(sub_feature) == 0) { return(NULL) }
                .stat_lipidTrend(
                    X=X[, sub_feature], X.info=subset(X.info, sub_feature),
                    group, radius, own_contri, test, permute_time, abun_weight)
            })
            names(output.df) <- names(chain_groups)
            return(new(
                "LipidTrendSE", lipid_se, split_chain=TRUE,
                even_chain_result=output.df$even,
                odd_chain_result=output.df$odd))
        } else {
            stop("chain_col is not present in the column name of rowData.")
        }
    } else if (isFALSE(split_chain)) {
        output.df <- .stat_lipidTrend(
            X, X.info, group, radius, own_contri, test, permute_time,
            abun_weight)
        return(new(
            "LipidTrendSE", lipid_se, split_chain=FALSE, result=output.df))
    } else {
        stop("split_chain must be logical value.")
    }
}

#' @title Region Statistics Calculation
#' @description Calculate region statistics using either t-test or Wilcoxon test
#' @param X Matrix. Lipid abundance data matrix.
#' @param Y Matrix. Group information matrix.
#' @param test Character. Statistical test to use ("t.test" or "Wilcoxon").
#' Default is "t.test".
#' @return Numeric vector of statistical test results.
#' @keywords internal
.regionStat <- function(X, Y, test="t.test", ...){
    n1 <- sum(Y[,1])
    n2 <- sum(1-Y[,1])
    G1 <- Y / n1
    G2 <- (1-Y) / n2
    mean.group1 <- t(X) %*% G1
    mean.group2 <- t(X) %*% G2
    if (test == "t.test") {
        var.group1 <- (t(X^2) %*% Y  - n1*mean.group1^2)/(n1 - 1)
        var.group2 <- (t(X^2) %*% (1-Y)  - n2*mean.group2^2)/(n2 - 1)
        pooled.var <- ((n1-1)*var.group1 + (n2-1)*var.group2)/(n1+n2-2)
        t.stat <- (mean.group1-mean.group2) / sqrt((1/n1+1/n2)*pooled.var)
        df <- n1 + n2 - 2
        t.test.pval <- 2*pt(abs(t.stat), df = df, lower.tail = FALSE)
        res <- -log10(t.test.pval) * sign(t.stat)
    }
    if (test == "Wilcoxon") {
        wilcox.pval <- apply(
            Y, 2, function(y){apply(
                X, 2, function(x) wilcox.test(
                    x ~ y, alternative="two.sided")$p.value)
            })
        res <- -log10(wilcox.pval) * sign(mean.group1-mean.group2)
    }
    return(res)
}

#' @title Statistical Analysis for lipidTrend
#' @description Perform statistical analysis for lipid trend data
#' @param X Matrix. Lipid abundance data matrix.
#' @param X.info Data frame. Lipid feature information.
#' @param group Vector. Group assignments.
#' @param radius Numeric. Distance of neighbor to include.
#' @param own_contri Numeric. Proportion of own contribution.
#' @param test Character. Statistical test to use.
#' @param permute_time Integer. Number of permutations.
#' @param abun_weight Logical. Whether to use abundance weighting.
#' @return Data frame containing analysis results.
#' @keywords internal
.stat_lipidTrend <- function(
        X, X.info, group, radius, own_contri, test, permute_time, abun_weight){
    region.stat.obs <- .regionStat(X=X, Y=as.matrix(group), test=test)
    dis_res <- .countDistance(X.info)
    dist_input <- dis_res$normalized_matrix
    x.distance <- dis_res$x_distance
    y.distance <- dis_res$y_distance
    dim <- dis_res$dimensions
    dist.mat <- as.matrix(dist(dist_input, method = "euclidean"))
    dist.mat[dist.mat>=radius] <- Inf
    stat_res <- .smooth_permutation(
        X, group, dist.mat, own_contri, test, abun_weight, permute_time,
        region.stat.obs)
    smooth.stat <- stat_res$smooth.stat
    smooth.stat.permute <- stat_res$smooth.stat.permute
    output.df <- .resultSummary(
        X, X.info, group, smooth.stat, smooth.stat.permute, region.stat.obs,
        dim, permute_time)
    return(output.df)
}

#' @title Smooth Permutation Analysis
#' @description Perform smooth permutation analysis for lipid trend data
#' @param X Matrix. Lipid abundance data matrix.
#' @param group Vector. Group assignments.
#' @param dist.mat Matrix. Distance matrix.
#' @param own_contri Numeric. Proportion of own contribution.
#' @param test Character. Statistical test to use.
#' @param abun_weight Logical. Whether to use abundance weighting.
#' @param permute_time Integer. Number of permutations.
#' @param region.stat.obs Numeric vector. Observed region statistics.
#' @return List containing smooth statistics and permutation results.
#' @keywords internal
.smooth_permutation <- function(
        X, group, dist.mat, own_contri, test, abun_weight, permute_time,
        region.stat.obs) {
    dist.idx <- which.max(colSums(exp(-dist.mat^2)))
    radius.seq <- seq(0.01, 5, 0.01)
    own_contri.seq <- vapply(
        radius.seq, function(x) 1/sum(exp(-dist.mat[,dist.idx]^2*x)),
        FUN.VALUE=numeric(1))
    radius.idx <- sum(own_contri.seq <= own_contri)
    if (radius.idx + 1 > length(radius.seq)) {
        stop(sprintf(
            "Based on the data structure and 'radius' setting, \n
        the 'own_contri' should smaller than %0.4f", max(own_contri.seq)))
    }
    kernel.radius <- radius.seq[radius.idx + 1]
    exp.dist <- exp(-dist.mat^2*kernel.radius)
    rescaled.X <- X-min(X)
    rescaled.X <- 4*rescaled.X/max(rescaled.X) + 1
    avg.log.expr <- colMeans(rescaled.X)
    if (abun_weight) {
        smooth.stat <- exp.dist %*% (region.stat.obs * avg.log.expr)
    } else {
        smooth.stat <- exp.dist %*% region.stat.obs
    }
    Y.permute <- vapply(
        seq_len(permute_time), function(x) sample(group, replace=FALSE),
        FUN.VALUE=numeric(length(group)))
    region.stat.permute <- .regionStat(X=X, Y=Y.permute, test=test)
    if (abun_weight) {
        smooth.stat.permute <- exp.dist %*% (region.stat.permute * avg.log.expr)
    } else {
        smooth.stat.permute <- exp.dist %*% region.stat.permute
    }
    return(list(
        smooth.stat=smooth.stat, smooth.stat.permute=smooth.stat.permute))
}

#' @title Summarize Results
#' @description Summarize and format lipid trend analysis results
#' @param X Matrix. Lipid abundance data matrix.
#' @param X.info Data frame. Lipid feature information.
#' @param group Vector. Group assignments.
#' @param smooth.stat Numeric vector. Smoothed statistics.
#' @param smooth.stat.permute Matrix. Permutation results.
#' @param region.stat.obs Numeric vector. Observed region statistics.
#' @param dimension Integer. Dimension of the analysis.
#' @param permute_time Integer. Number of permutations performed.
#' @param ... Additional arguments.
#' @return Data frame containing formatted results.
#' @importFrom dplyr mutate
#' @importFrom stats p.adjust
#' @keywords internal
.resultSummary <- function(
        X, X.info, group, smooth.stat, smooth.stat.permute,
        region.stat.obs, dimension, permute_time, ...){
    smooth.stat.all <- cbind(smooth.stat, smooth.stat.permute)
    smooth.stat.pval <- apply(
        smooth.stat.all, 1, function(x) mean(abs(x[1]) < abs(x[-1])))
    smooth.stat.pval[smooth.stat.pval == 0] <- 1/permute_time
    smooth.stat.pval.BH <- p.adjust(smooth.stat.pval, method="BH")
    direction <- ifelse(sign(smooth.stat) > 0, "+", "-")
    marginal.pval <- 10^-abs(region.stat.obs)
    marginal.pval.BH <- p.adjust(marginal.pval, method="BH")
    log2.FC <- apply(
        X, 2, function(x) log2(mean(x[group == 1])/mean(x[group == 0])))
    FC.direction <- ifelse(log2.FC > 0, "+", "-")
    smooth.stat.pval.BH[FC.direction != direction] <- 1
    result <- data.frame(
        X.info, direction=direction, smoothing.pval.BH=smooth.stat.pval.BH,
        marginal.pval.BH=marginal.pval.BH, log2.FC=log2.FC)
    if (dimension == 2) {
        result %<>%
            dplyr::mutate(avg.expr=round(colMeans(X), 2), .before="direction")
    }
    if (dimension == 1) {
        result %<>%
            dplyr::mutate(
                avg.expr.ctrl=round(colMeans(X[group == 0,]), 2),
                avg.expr.case=round(colMeans(X[group == 1,]), 2),
                .before="direction")
    }
    return(result)
}
