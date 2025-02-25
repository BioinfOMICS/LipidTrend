library(testthat)
library(SummarizedExperiment)

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

test_that("plotRegion1D creates valid plot", {
    set.seed(1234)
    se <- create_mock_se()
    se_split <- create_mock_se(split_chain=TRUE)
    result <- analyzeLipidRegion(
        lipid_se=se, ref_group="control", split_chain=FALSE, chain_col=NULL,
        radius=3,own_contri=0.5, permute_time=100)
    expect_no_error(plot <- plotRegion1D(result, p_cutoff=0.05))
    expect_true(inherits(plot, "ggplot"))
    # positive regions
    result_df <- getResult(result)
    result_df$direction <- "+"
    result_df$smoothing.pval.BH <- 0.01
    result_df$avg.expr.ctrl <- 1
    result_df$avg.expr.case <- 2
    attr(result, "result") <- result_df
    expect_no_error(plot_pos <- plotRegion1D(result, p_cutoff=0.05))
    expect_true(inherits(plot_pos, "ggplot"))
    # negative regions
    result_df$direction <- "-"
    result_df$avg.expr.ctrl <- 2
    result_df$avg.expr.case <- 1
    attr(result, "result") <- result_df
    expect_no_error(plot_neg <- plotRegion1D(result, p_cutoff=0.05))
    expect_true(inherits(plot_neg, "ggplot"))
    # mixed positive and negative regions
    result_df$direction <- rep(c("+", "-"), length.out=nrow(result_df))
    result_df$smoothing.pval.BH <- rep(c(0.01, 0.1), length.out=nrow(result_df))
    attr(result, "result") <- result_df
    expect_no_error(plot_mixed <- plotRegion1D(result, p_cutoff=0.05))
    expect_true(inherits(plot_mixed, "ggplot"))
})

test_that("plotRegion1D handles split chain analysis correctly", {
    set.seed(1234)
    se_split <- create_mock_se(split_chain=TRUE)
    result <- analyzeLipidRegion(
        lipid_se=se_split, ref_group="control", split_chain=TRUE,
        chain_col="chain", radius=3, own_contri=0.5, permute_time=100)
    # both even and odd chain
    expect_no_error(plot <- plotRegion1D(result, p_cutoff=0.05))
    expect_true(inherits(plot$even_result, "ggplot"))
    expect_true(inherits(plot$odd_result, "ggplot"))
    # one chain type has no significant regions
    even_results <- getEvenChainResult(result)
    even_results$smoothing.pval.BH <- 1
    attr(result, "even_chain_results") <- even_results
    expect_no_error(plot_one_sig <- plotRegion1D(result, p_cutoff=0.05))
})

test_that("plotRegion1D handles invalid input and edge cases", {
    se <- create_mock_se()
    result <- analyzeLipidRegion(
        lipid_se=se, ref_group="control", split_chain=FALSE, chain_col=NULL,
        radius=3, own_contri=0.5, permute_time=100)
    # invalid p_cutoff values
    expect_error(
        plotRegion1D(result, p_cutoff=-1),
        "p_cutoff must be a numeric value between 0 and 1"
    )
    expect_error(
        plotRegion1D(result, p_cutoff=2),
        "p_cutoff must be a numeric value between 0 and 1"
    )
    # no significant regions
    result_df <- getResult(result)
    result_df$smoothing.pval.BH <- 1
    attr(result, "result") <- result_df
    expect_no_error(plot_no_sig <- plotRegion1D(result, p_cutoff=0.05))
    expect_true(inherits(plot_no_sig, "ggplot"))
    # missing results
    bad_result <- result
    attr(bad_result, "result") <- NULL
    expect_error(plotRegion1D(bad_result, p_cutoff=0.05))
    # invalid split chain status
    bad_split_result <- result
    attr(bad_split_result, "split_chain") <- "invalid"
    expect_error(plotRegion1D(bad_split_result, p_cutoff=0.05))
})

test_that("plotRegion1D handles various p-value cutoffs", {
    se <- create_mock_se()
    result <- analyzeLipidRegion(
        lipid_se=se, ref_group="control", split_chain=FALSE, permute_time=100)
    p_cutoffs <- c(0.01, 0.05, 0.1)
    for(p_cut in p_cutoffs) {
        expect_no_error(plot <- plotRegion1D(result, p_cutoff=p_cut))
        expect_true(inherits(plot, "ggplot"))
    }
})
