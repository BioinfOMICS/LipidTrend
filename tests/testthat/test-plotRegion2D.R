library(testthat)
library(SummarizedExperiment)

# Helper function to create mock SummarizedExperiment object
create_mock_se <- function(n_samples=6, n_features=10, split_chain=FALSE) {
    if (split_chain) {
        feature_names <- as.character(seq(33, 33 + n_features - 1, by = 1))
        assay_data <- matrix(
            rnorm(n_samples * n_features, mean=5, sd=1),
            nrow=n_features, ncol=n_samples, dimnames=list(feature_names, NULL))
        row_data <- data.frame(
            chain=as.numeric(feature_names), x=as.numeric(feature_names),
            row.names=feature_names)
    } else {
        base_chains <- seq(33, 33 + n_features - 1)
        double_bonds <- rep(1:2, length.out = n_features)
        feature_names <- paste("TG", paste0(base_chains, ":", double_bonds))
        assay_data <- matrix(
            rnorm(n_samples * n_features, mean=5, sd=1),
            nrow=n_features, ncol=n_samples, dimnames=list(feature_names, NULL))
        row_data <- data.frame(
            x=base_chains, y=double_bonds, Total.C=base_chains,
            Total.DB=double_bonds, row.names=feature_names)
    }
    col_data <- data.frame(
        sample_name=paste0("sample", seq_len(n_samples)),
        label_name=paste0("label", seq_len(n_samples)),
        group=rep(c("control", "case"), length.out=n_samples),
        row.names=paste0("sample", seq_len(n_samples)))
    SummarizedExperiment(
        assays=list(abundance=assay_data), rowData=row_data, colData=col_data)
}

test_that("plotRegion2D creates valid plot", {
    set.seed(1234)
    n_features <- 25
    se <- create_mock_se(n_features=n_features)
    x_coords <- rep(seq_len(5), each=5)
    y_coords <- rep(seq_len(5), times=5)
    base_chains <- 33:(33 + n_features - 1)
    feature_names <- paste("TG", paste0(base_chains, ":", y_coords))
    row_data <- data.frame(
        Total.C=x_coords, Total.DB=y_coords, row.names=feature_names)
    rowData(se) <- row_data
    result <- analyzeLipidRegion(
        lipid_se=se, ref_group="control", split_chain=FALSE,
        chain_col=NULL, radius=3, own_contri=0.5, permute_time=100)
    expect_no_error(plot <- plotRegion2D(result, p_cutoff=0.05))
    expect_true(inherits(plot, "ggplot"))
    # high regions
    result_df <- getResult(result)
    result_df$direction <- "+"
    result_df$smoothing.pval.BH <- 0.01
    result_df$log2.FC <- 2
    attr(result, "result") <- result_df
    expect_no_error(plot_high <- plotRegion2D(result, p_cutoff=0.05))
    expect_true(inherits(plot_high, "ggplot"))
    # low regions
    result_df$direction <- "-"
    result_df$log2.FC <- -2
    attr(result, "result") <- result_df
    expect_no_error(plot_low <- plotRegion2D(result, p_cutoff=0.05))
    expect_true(inherits(plot_low, "ggplot"))
})

test_that("plotRegion2D handles split chain analysis correctly", {
    set.seed(1234)
    n_features <- 25
    se_split <- create_mock_se(n_features=n_features, split_chain=TRUE)
    x_coords <- rep(seq_len(5), each=5)
    y_coords <- rep(seq_len(5), times=5)
    row_data <- data.frame(
        Total.C=x_coords, Total.DB=y_coords,
        row.names=rownames(assay(se_split)))
    rowData(se_split) <- row_data
    result <- analyzeLipidRegion(
        lipid_se=se_split, ref_group="control", split_chain=TRUE,
        chain_col="Total.C", radius=3, own_contri=0.5, permute_time=100)
    # split chain
    expect_no_error(plot <- plotRegion2D(result, p_cutoff=0.05))
    expect_true(inherits(plot$even_result, "ggplot"))
    expect_true(inherits(plot$odd_result, "ggplot"))
    # one chain type has no significant regions
    even_results <- getEvenChainResult(result)
    even_results$smoothing.pval.BH <- 1
    attr(result, "even_chain_results") <- even_results
    expect_no_error(plot_partial <- plotRegion2D(result, p_cutoff=0.05))
})

test_that("plotRegion2D handles wall building correctly", {
    set.seed(1234)
    n_features <- 25
    se <- create_mock_se(n_features=n_features)
    x_coords <- rep(seq_len(5), each=5)
    y_coords <- rep(seq_len(5), times=5)
    row_data <- data.frame(
        Total.C=x_coords, Total.DB=y_coords, row.names=rownames(assay(se)))
    rowData(se) <- row_data
    result <- analyzeLipidRegion(
        lipid_se=se, ref_group="control", split_chain=FALSE, radius=3,
        own_contri=0.5, permute_time=100)
    result_df <- getResult(result)
    result_df$direction <- rep(c("+", "-"), length.out=nrow(result_df))
    result_df$smoothing.pval.BH <- 0.01
    result_df$log2.FC <- ifelse(result_df$direction == "+", 2, -2)
    attr(result, "result") <- result_df
    expect_no_error(plot <- plotRegion2D(result, p_cutoff=0.05))
    expect_true(inherits(plot, "ggplot"))
    # with no walls (all non-significant)
    result_df$smoothing.pval.BH <- 1
    attr(result, "result") <- result_df
    expect_no_error(plot_no_walls <- plotRegion2D(result, p_cutoff=0.05))
    expect_true(inherits(plot_no_walls, "ggplot"))
    # with all walls (all significant with alternating regions)
    result_df$smoothing.pval.BH <- 0.01
    result_df$direction <- rep(c("+", "-"), length.out=nrow(result_df))
    attr(result, "result") <- result_df
    expect_no_error(plot_all_walls <- plotRegion2D(result, p_cutoff=0.05))
    expect_true(inherits(plot_all_walls, "ggplot"))
    # with specific region patterns
    result_df$direction[c(1:5, 21:25)] <- "+"
    result_df$direction[6:20] <- "-"
    result_df$smoothing.pval.BH <- 0.01
    attr(result, "result") <- result_df
    expect_no_error(plot_patterns <- plotRegion2D(result, p_cutoff=0.05))
    expect_true(inherits(plot_patterns, "ggplot"))
    # with different p-value annotations
    result_df$marginal.pval.BH <- rep(c(0.0001, 0.001, 0.01, 0.05, 0.1), 5)
    attr(result, "result") <- result_df
    expect_no_error(plot_pvals <- plotRegion2D(result, p_cutoff=0.05))
    expect_true(inherits(plot_pvals, "ggplot"))
    # with varied abundance values
    result_df$avg.expr <- seq(1, 25)
    attr(result, "result") <- result_df
    expect_no_error(plot_abundance <- plotRegion2D(result, p_cutoff=0.05))
    expect_true(inherits(plot_abundance, "ggplot"))
})

test_that("plotRegion2D handles invalid inputs and edge cases", {
    se <- create_mock_se()
    result <- analyzeLipidRegion(
        lipid_se=se, ref_group="control", split_chain=FALSE, radius=3,
        own_contri=0.5, permute_time=100)
    # invalid p_cutoff
    expect_error(
        plotRegion2D(result, p_cutoff=2),
        "p_cutoff must be a numeric value between 0 and 1"
    )
    expect_error(
        plotRegion2D(result, p_cutoff=-1),
        "p_cutoff must be a numeric value between 0 and 1"
    )
    # invalid log2FC_cutoff
    expect_error(
        plotRegion2D(result, p_cutoff=0.05, log2FC_cutoff=-1),
        "log2FC_cutoff must be a positive numeric value"
    )
    # invalid split chain status
    bad_result <- result
    attr(bad_result, "split_chain") <- "invalid"
    expect_error(plotRegion2D(bad_result, p_cutoff=0.05))
    # missing required columns
    result_df <- getResult(result)
    result_df$Total.C <- NULL
    attr(result, "result") <- result_df
    expect_error(plotRegion2D(result, p_cutoff=0.05))
})

test_that("plotRegion2D handles irregular grids", {
    n_features <- 15
    se <- create_mock_se(n_features=n_features)
    x_coords <- rep(seq_len(5), each=3)
    y_coords <- rep(seq_len(3), times=5)
    row_data <- data.frame(
        Total.C=x_coords, Total.DB=y_coords, row.names=rownames(assay(se)))
    rowData(se) <- row_data
    result <- analyzeLipidRegion(
        lipid_se=se, ref_group="control", split_chain=FALSE, radius=3,
        own_contri=0.5, permute_time=100)
    expect_no_error(plot <- plotRegion2D(result, p_cutoff=0.05))
    expect_true(inherits(plot, "ggplot"))
})

test_that(".buildWall validates input parameters", {
    # Test non-numeric feature.idx
    expect_error(
        .buildWall(
            feature.idx="1", X.info=data.frame(x=seq_len(3), y=seq_len(3)),
            selected.region=c("High", "Low", "None"),
            x.distance=1,y.distance=1),
        "feature.idx must be a single numeric value"
    )
    # Test length > 1 feature.idx
    expect_error(
        .buildWall(
            feature.idx=c(1, 2), X.info=data.frame(x=seq_len(3), y=seq_len(3)),
            selected.region=c("High", "Low", "None"),
            x.distance=1, y.distance=1),
        "feature.idx must be a single numeric value"
    )
    # Test non-numeric distance parameters
    expect_error(
        .buildWall(
            feature.idx=1, X.info=data.frame(x=seq_len(3), y=seq_len(3)),
            selected.region=c("High", "Low", "None"),
            x.distance="1", y.distance=1),
        "distance parameters must be numeric"
    )
    expect_error(
        .buildWall(
            feature.idx=1, X.info=data.frame(x=seq_len(3), y=seq_len(3)),
            selected.region=c("High", "Low", "None"),
            x.distance=1, y.distance="1"),
        "distance parameters must be numeric"
    )
})

test_that(".pvalAnnotation handles edge cases", {
    # Test empty input
    expect_equal(.pvalAnnotation(numeric(0)), "")
    # Test NA input
    expect_equal(.pvalAnnotation(NA), "")
    # Test all p-value ranges
    expect_equal(.pvalAnnotation(0.0001), "***")
    expect_equal(.pvalAnnotation(0.001), "**")
    expect_equal(.pvalAnnotation(0.009), "**")
    expect_equal(.pvalAnnotation(0.04), "*")
    expect_equal(.pvalAnnotation(0.09), ".")
    expect_equal(.pvalAnnotation(0.5), "")
})

test_that(".plotHeatmap validates input and creates plot", {
    # Create test data
    X.info <- data.frame(x=rep(seq_len(3), each=3), y=rep(seq_len(3), times=3))
    # Test missing required columns
    heatmap.df <- data.frame(v1=seq_len(9), v2=seq_len(9))
    expect_error(
        .plotHeatmap(X.info, heatmap.df, NULL, 1, 1),
        "Required columns missing from heatmap.df"
    )
    # Test invalid distance parameters
    heatmap.df <- data.frame(
        v1=seq_len(9), v2=seq_len(9), avg.expr=runif(9), log2.FC=rnorm(9),
        pval.annotate=rep("*", 9))
    expect_error(
        .plotHeatmap(X.info, heatmap.df, NULL, "1", 1),
        "Distance parameters must be numeric"
    )
})

test_that(".calcSegment handles all wall directions", {
    # Create test wall data
    wall_data <- data.frame(x=2, y=2)
    x_dist <- 1
    y_dist <- 1
    # Test top wall
    top_segment <- .calcSegment(wall_data, "top", x_dist, y_dist)
    expect_equal(
        top_segment,
        list(
            x=wall_data$x - x_dist/2, xend=wall_data$x + x_dist/2,
            y=wall_data$y + y_dist/2, yend=wall_data$y + y_dist/2)
    )
    # Test right wall
    right_segment <- .calcSegment(wall_data, "right", x_dist, y_dist)
    expect_equal(
        right_segment,
        list(
            x=wall_data$x + x_dist/2, xend=wall_data$x + x_dist/2,
            y=wall_data$y - y_dist/2, yend=wall_data$y + y_dist/2)
    )
    # Test bottom wall
    bottom_segment <- .calcSegment(wall_data, "bottom", x_dist, y_dist)
    expect_equal(
        bottom_segment,
        list(
            x=wall_data$x - x_dist/2, xend=wall_data$x + x_dist/2,
            y=wall_data$y - y_dist/2, yend=wall_data$y - y_dist/2)
    )
    # Test left wall
    left_segment <- .calcSegment(wall_data, "left", x_dist, y_dist)
    expect_equal(
        left_segment,
        list(
            x=wall_data$x - x_dist/2, xend=wall_data$x - x_dist/2,
            y=wall_data$y - y_dist/2,yend=wall_data$y + y_dist/2)
    )
    # Test with NULL wall data
    expect_null(.calcSegment(NULL, "top", x_dist, y_dist))
})
