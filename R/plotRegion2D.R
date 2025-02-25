#' @title Plot region trends for two-dimensional lipid features
#' @description plotRegion2D visualizes lipid trend analysis results in two
#' dimensions. The returned heatmap will show the results.
#' @param object A LipidTrendSE object containing analysis results.
#' @param p_cutoff Numeric. Significance level. The range is 0 to 1.
#'   Default is \code{0.05}.
#' @param log2FC_cutoff Numeric. Significant level of log2-transformed fold
#' change. Default is \code{3}.
#' @return For split chain analyses, a list with two elements:
#'   \item{even_result}{Heatmap for even chain features, or NULL if none exist}
#'   \item{odd_result}{Heatmap for odd chain features, or NULL if none exist}
#'   For non-split analyses, a single heatmap.
#' @examples
#' data("lipid_se_2D")
#' res_se <- analyzeLipidRegion(
#'     lipid_se=lipid_se_2D, ref_group="sgCtrl", split_chain=FALSE,
#'     chain_col=NULL, radius=3, own_contri=0.5, permute_time=100,
#'     abun_weight=TRUE)
#' plot_2D <- plotRegion2D(res_se, p_cutoff=0.05)
#' @importFrom ggplot2 ggplot aes geom_point theme_bw scale_colour_gradient2
#'   geom_text scale_x_continuous scale_y_continuous labs scale_size
#'   geom_segment
#' @export
#' @seealso
#' \code{\link{analyzeLipidRegion}} for generating the input LipidTrendSE object
setGeneric("plotRegion2D", function(object, p_cutoff=0.05, log2FC_cutoff=3)
    standardGeneric("plotRegion2D"))

#' @rdname plotRegion2D
#' @aliases plotRegion2D,LipidTrendSE-method
#' @exportMethod plotRegion2D
setMethod(
    "plotRegion2D", "LipidTrendSE",
    function(object, p_cutoff=0.05, log2FC_cutoff=3) {
    if (!is.numeric(p_cutoff) || p_cutoff < 0 || p_cutoff > 1) {
        stop("p_cutoff must be a numeric value between 0 and 1")
    }
    if (!is.numeric(log2FC_cutoff) || log2FC_cutoff <= 0) {
        stop("log2FC_cutoff must be a positive numeric value")
    }
    if (.getSplitChainStatus(object)) {
        even_results <- getEvenChainResult(object)
        odd_results <- getOddChainResult(object)
        plots <- list(
            even_result=if(
                !is.null(even_results)) .trendPlot2D(
                    even_results, p_cutoff, log2FC_cutoff) else NULL,
            odd_result=if(
                !is.null(odd_results)) .trendPlot2D(
                    odd_results, p_cutoff, log2FC_cutoff) else NULL
        )
        return(plot2D=plots)
    } else {
        results <- getResult(object)
        plot <- if(!is.null(results)) .trendPlot2D(
            results, p_cutoff, log2FC_cutoff) else NULL
        return(plot2D=plot)
    }
})

#' @title Create 2D Trend Plot
#' @description Creates a 2D visualization of lipid trends with
#' statistical annotations
#' @param res Data frame containing analysis results
#' @param p_cutoff Numeric significance level cutoff
#' @param log2FC_cutoff Numeric log2 fold change cutoff
#' @return A ggplot2 object representing the trend visualization
#' @keywords internal
.trendPlot2D <- function(res, p_cutoff, log2FC_cutoff){
    X.info <- res[, seq_len(2)]
    direction.int <- ifelse(res$direction == "+", 1, -1)
    smoothing.pval <- res$smoothing.pval.BH
    selected.regions <- list(
        high=ifelse(smoothing.pval<p_cutoff&direction.int>0, "High", "None"),
        low=ifelse(smoothing.pval<p_cutoff&direction.int<0, "Low", "None"))
    dis_res <- .countDistance(X.info)
    x.distance <- dis_res$x_distance
    y.distance <- dis_res$y_distance
    buildAllWalls <- function(region_type) {
        walls <- lapply(seq_len(nrow(X.info)), function(x) {
            .buildWall(
                feature.idx=x, X.info, selected.regions[[region_type]],
                x.distance, y.distance)
        })
        walls_df <- as.data.frame(do.call(rbind, walls))
        colnames(walls_df)[seq_len(2)] <- c("x","y")
        split(walls_df[, seq_len(2)], walls_df$pos)
    }
    walls <- list(high=buildAllWalls("high"),low=buildAllWalls("low"))
    heatmap.df <- data.frame(
        v1=X.info[,1], v2=X.info[,2], avg.expr=res$avg.expr,
        log2.FC=pmin(pmax(res$log2.FC, -log2FC_cutoff), log2FC_cutoff),
        pval.annotate=vapply(
            res$marginal.pval.BH, function(pval) .pvalAnnotation(pval),
            FUN.VALUE=character(1)),
        signed.logp.smooth=direction.int*smoothing.pval)
    heatmap <- .plotHeatmap(X.info, heatmap.df, walls, x.distance, y.distance)
    return(heatmap)
}

#' @title Build Walls Between Different Regions
#' @description Identifies walls needed between different regions based on
#'   spatial relationships and region assignments
#' @param feature.idx Numeric. Index of the feature to analyze
#' @param X.info Matrix or data.frame. Contains spatial coordinates
#' @param selected.region Character vector. Region assignments for features
#' @param x.distance Numeric. Scaling factor for x-coordinates
#' @param y.distance Numeric. Scaling factor for y-coordinates
#' @return Data frame with coordinates and wall positions
#' @importFrom stats setNames
#' @keywords internal
.buildWall <- function(
        feature.idx, X.info, selected.region, x.distance, y.distance) {
    if (!is.numeric(feature.idx) || length(feature.idx) != 1) {
        stop("feature.idx must be a single numeric value")
    }
    if (!is.numeric(x.distance) || !is.numeric(y.distance)) {
        stop("distance parameters must be numeric")
    }
    dis_res <- .countDistance(X.info)
    dist_matrix <- dis_res$normalized_matrix
    ref_point <- dist_matrix[feature.idx, ]
    neighbors <- lapply(
        c("right", "left", "top", "bottom"), function(direction) {
            .findNeighbor(ref_point, dist_matrix, direction)})
    names(neighbors) <- c("right", "left", "top", "bottom")
    wall_positions <- character()
    needsWall <- function(neighbor_indices, direction) {
        feature_region <- selected.region[feature.idx]
        n_neighbors <- sum(neighbor_indices)
        if (direction %in% c("left", "bottom")) {
            return(feature_region != "None" && n_neighbors == 0)
        }
        if (n_neighbors == 1) {
            neighbor_region <- selected.region[which(neighbor_indices)]
            return(feature_region != neighbor_region)
        } else {
            return(feature_region != "None")
        }
    }
    # for (direction in names(neighbors)) {
    #     if (needsWall(neighbors[[direction]], direction)) {
    #         wall_positions <- c(wall_positions, direction)
    #     }
    # }
    need_walls <- vapply(names(neighbors), function(direction) {
        needsWall(neighbors[[direction]], direction)
    }, logical(1))
    wall_positions <- names(neighbors)[need_walls]
    if (length(wall_positions) == 0) {
        wall_out <- data.frame(matrix(ncol = ncol(X.info) + 1, nrow = 0))
    } else {
        wall_out <- data.frame(do.call(
            rbind, replicate(
                length(wall_positions), as.numeric(X.info[feature.idx, ]),
                simplify = FALSE)), pos=wall_positions)
    }
    colnames(wall_out) <- c(colnames(X.info), "pos")
    return(wall_out)
}

#' @title Find Neighboring Features
#' @description Identifies neighboring features in specified directions
#' @param self Numeric vector. Reference point coordinates
#' @param others Matrix. Coordinates of other features
#' @param type Character. Direction to check ("right", "left", "top", "bottom")
#' @return Logical vector indicating neighbor presence
#' @keywords internal
.findNeighbor <- function(
        self, others, type=c("right", "left", "top", "bottom")) {
    type <- match.arg(type)
    x_diff <- others[, 1] - self[1]
    y_diff <- others[, 2] - self[2]
    switch(
        type,
        "right"  = y_diff == 0 & x_diff == 1,
        "left"   = y_diff == 0 & x_diff == -1,
        "top"    = x_diff == 0 & y_diff == 1,
        "bottom" = x_diff == 0 & y_diff == -1)
}

#' @title Create P-value Annotation
#' @description Converts p-values to significance symbols
#' @param pval Numeric. P-value to convert
#' @return Character string with significance symbols
#' @keywords internal
.pvalAnnotation <- function(pval) {
    if (length(pval) == 0 || is.na(pval)) return("")
    dplyr::case_when(
        pval < 0.001 ~ "***", pval < 0.01 ~ "**", pval < 0.05 ~ "*",
        pval < 0.1 ~ ".", TRUE ~ "")
}

#' @title Plot Heatmap with Abundance Data and Region Walls
#' @description Creates a ggplot2 heatmap visualization showing abundance data
#'   with walls between regions, using vectorized operations.
#' @param X.info Matrix or data.frame. Contains coordinate information
#' @param heatmap.df Data frame. Contains abundance and fold change data
#' @param walls List. Contains wall coordinates for different regions
#' @param x.distance Numeric. Distance between x-axis points
#' @param y.distance Numeric. Distance between y-axis points
#' @return A ggplot2 object
#' @importFrom ggplot2 ggplot aes geom_point theme_bw scale_colour_gradient2
#'   geom_text scale_x_continuous scale_y_continuous labs scale_size
#'   geom_segment
#' @keywords internal
.plotHeatmap <- function(X.info, heatmap.df, walls, x.distance, y.distance) {
    if (
        !all(
            c("v1", "v2", "avg.expr", "log2.FC", "pval.annotate") %in% colnames(
                heatmap.df))) {
        stop("Required columns missing from heatmap.df")
    }
    if (!is.numeric(x.distance) || !is.numeric(y.distance)) {
        stop("Distance parameters must be numeric")
    }
    base_plot <- ggplot(heatmap.df,aes(x=v1,y=v2)) + theme_bw() +
        geom_point(aes(size=avg.expr, color=log2.FC)) +
        scale_colour_gradient2(
            high="#F8766D", mid="white", low="blue", midpoint=0) +
        geom_text(aes(label=pval.annotate), color="black", size=4) +
        scale_x_continuous(
            breaks=seq(min(X.info[,1]), max(X.info[,1]), x.distance)) +
        scale_y_continuous(
            breaks=seq(min(X.info[,2]), max(X.info[,2]), y.distance)) +
        labs(x=colnames(X.info)[1], y=colnames(X.info)[2]) +
        scale_size(range=c(2,10))
    regions <- c("high", "low")
    directions <- c("top", "right", "bottom", "left")
    colors <- c(high = "red", low = "blue")
    segment_layers <- unlist(
        lapply(regions, function(region) {
            vapply(directions, function(direction) {
                layer <- .createSegmentLayer(
                    walls[[region]][[direction]], direction, colors[region],
                    x.distance, y.distance)
                if (is.null(layer)) {
                    return(list(NULL))
                } else {
                    return(list(layer))
                }
            }, FUN.VALUE=list(NULL))
        }), recursive = FALSE)
    segment_layers <- Filter(Negate(is.null), segment_layers)
    if (length(segment_layers) > 0) {
        final_plot <- Reduce(`+`, segment_layers, init=base_plot)
    } else {
        final_plot <- base_plot
    }
    return(final_plot)
}

#' @title Create Segment Layer for Wall Visualization
#' @description Creates a ggplot2 layer for wall segments
#' @param wall_data Data frame. Wall coordinates
#' @param direction Character. Wall direction
#' @param color Character. Wall color
#' @param x_dist Numeric. X-axis distance scale
#' @param y_dist Numeric. Y-axis distance scale
#' @return A ggplot2 layer object or NULL
#' @keywords internal
.createSegmentLayer <- function(wall_data, direction, color, x_dist, y_dist){
    if (is.null(wall_data)) {return(NULL)}
    segments <- .calcSegment(wall_data, direction, x_dist, y_dist)
    segment_df <- data.frame(
        x=segments$x, xend=segments$xend, y=segments$y, yend=segments$yend)
    geom_segment(
        data=segment_df, mapping=aes(x=x, xend=xend, y=y, yend=yend),
        linewidth=1, colour=color)
}

#' @title Calculate Segment Coordinates
#' @description Calculates start and end coordinates for wall segments
#' @param wall_data Data frame. Wall position data
#' @param direction Character. Wall direction
#' @param x_dist Numeric. X-axis distance scale
#' @param y_dist Numeric. Y-axis distance scale
#' @return List of segment coordinates
#' @keywords internal
.calcSegment <- function(wall_data, direction, x_dist, y_dist) {
    if (is.null(wall_data)) {return(NULL)}
    switch(
        direction,
        "top"=list(
            x=wall_data$x-x_dist/2, xend=wall_data$x+x_dist/2,
            y=wall_data$y+y_dist/2, yend=wall_data$y+y_dist/2),
        "right"=list(
            x=wall_data$x+x_dist/2, xend=wall_data$x+x_dist/2,
            y=wall_data$y-y_dist/2, yend=wall_data$y+y_dist/2),
        "bottom"=list(
            x=wall_data$x-x_dist/2, xend=wall_data$x+x_dist/2,
            y=wall_data$y-y_dist/2, yend=wall_data$y-y_dist/2),
        "left"=list(
            x=wall_data$x-x_dist/2, xend=wall_data$x-x_dist/2,
            y=wall_data$y-y_dist/2, yend=wall_data$y+y_dist/2))
}
