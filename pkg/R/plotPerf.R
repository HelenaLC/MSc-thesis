#' @rdname plotPerf
#' @title TPR vs. FDR curve
#'
#' @description 
#' Plots a TPR vs. FDR curve from the performance output of \code{iCOBRA}.
#'
#' @param x output from \code{iCOBRA::calculate_performance}
#' @param main character string to use as plot title.
#' @param cols named colors palette to use for methods. 
#' Defaults to \code{RColorBrewer::brewer.pal(11, "Spectral")},
#'
#' @author Helena Lucia Crowell
#'
#' @include utils.R
#' @import ggplot2
#' @importFrom grDevices colorRampPalette
#' @importFrom iCOBRA COBRAPerformance
#' @importFrom RColorBrewer brewer.pal
#'
#' @export

setMethod(f = "plotPerf", 
    signature = signature(x = "COBRAPerformance"), 
    definition = function(x, main = NULL, cols = NULL) {
        df <- data.frame(x@fdrtprcurve)
        methods <- unique(df$method)
        
        # get method colors
        n_methods <- length(methods)
        if (is.null(cols)) {
            cols <- brewer.pal(11, "Spectral")
            if (n_methods > 11)
                cols <- colorRampPalette(cols)(n_methods)
            cols <- cols[sample(n_methods)]
            cols <- setNames(cols, methods)
        }
        
        df$method <- factor(df$method, levels=methods)
        df <- df[!(is.na(df$FDR) | is.na(df$TPR)), ]
        ths <- c(.01,.05,.1)
        ggplot(df, aes_string(x="FDR", y="TPR", col="method")) +
            scale_color_manual(NULL, values=cols[methods]) +
            guides(color = guide_legend(NULL)) +
            geom_vline(xintercept=seq(0,1,.1), linetype=2, color="lightgrey", size=.2) +
            geom_vline(xintercept=ths, linetype=2, size=.4) +
            scale_x_continuous(limits=c(0,1), breaks=seq(0,1,.2), expand=c(.05,0)) +
            scale_y_continuous(limits=c(0,1), breaks=seq(0,1,.2), expand=c(.05,0)) +
            geom_path(size=.5) + ggtitle(main) + theme_bw() + prettify
    }
)