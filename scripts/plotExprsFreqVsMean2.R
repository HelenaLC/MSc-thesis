plotExprsFreqVsMean2 <- function (object, freq_exprs, mean_exprs, controls, by_show_single = FALSE, 
    show_smooth = TRUE, show_se = TRUE, ...) {
    if (!is(object, "SingleCellExperiment")) {
        stop("Object must be an SingleCellExperiment")
    }
    if (missing(freq_exprs) || is.null(freq_exprs)) {
        freq_exprs <- scater:::.qc_hunter(object, "n_cells_by_counts", 
            mode = "row")
    }
    if (missing(mean_exprs) || is.null(mean_exprs)) {
        mean_exprs <- scater:::.qc_hunter(object, "mean_counts", mode = "row")
    }
    if (missing(controls)) {
        controls <- scater:::.qc_hunter(object, "is_feature_control", 
            mode = "row")
    }
    freqs <- scater:::.choose_vis_values(object, freq_exprs, mode = "row", 
        search = "metadata")$val/ncol(object) * 100
    means <- log2(scater:::.choose_vis_values(object, mean_exprs, mode = "row", 
        search = "metadata")$val + 1)
    plot_out <- plotRowData(object, y = data.frame(`Percentage of expressing cells` = freqs, 
        check.names = FALSE), x = data.frame(`Mean expression` = means, 
            check.names = FALSE), colour_by = controls, shape_by = controls, 
        ...) + scale_shape_manual(values = c(1, 17))
    if (!is.null(controls)) {
        conts <- scater:::.choose_vis_values(object, controls, mode = "row", 
            search = "metadata", discard_solo = !by_show_single)$val
    } else {
        conts <- logical(nrow(object))
    }
    mn_vs_fq <- data.frame(mn = means, fq = freqs/100)
    text_x_loc <- min(mn_vs_fq$mn) + 0.6 * diff(range(mn_vs_fq$mn))
    if (show_smooth) {
        if (any(conts)) {
            tmp_tab <- mn_vs_fq[conts, ]
        } else {
            tmp_tab <- mn_vs_fq
        }
        plot_out <- plot_out + geom_smooth(aes_string(x = "mn", 
            y = "100 * fq"), data = tmp_tab, colour = "firebrick", 
            size = 1, se = show_se)
    }
    if (any(conts)) {
        dropout <- nls(fq ~ (1/(1 + exp(-(-i + 1 * mn)))), start = list(i = 5), 
            data = mn_vs_fq[conts, ])
        plot_out <- plot_out + geom_vline(xintercept = coef(dropout), 
            linetype = 2) + annotate("text", x = text_x_loc, 
                y = 60, size=2, label = paste(sum(mn_vs_fq[!conts, "mn"] > 
                        coef(dropout)), " genes with mean expression\nhigher than value for 50% dropout of feature controls", 
                    sep = ""))
    }
    plot_out <- plot_out + geom_hline(yintercept = 50, linetype = 2) + 
        annotate("text", x = text_x_loc, y = 40, size=2, label = paste(sum(mn_vs_fq$fq >= 
                0.5), " genes are expressed\nin at least 50% of cells", 
            sep = "")) + annotate("text", x = text_x_loc, y = 20, size=2,
                label = paste(sum(mn_vs_fq$fq >= 0.25), " genes are expressed\nin at least 25% of cells", 
                    sep = ""))
    plot_out
}
