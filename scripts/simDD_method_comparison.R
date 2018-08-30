# scDD, method comparison
# ==============================================================================
# Runs and compares (TPR vs. FDR curves) method performances 
# on a set of 12 simulated data sets, each containing only genes 
# from a single category, with increasing difficulty.
# ------------------------------------------------------------------------------

path <- "results/method_comparison"

# load packages
suppressPackageStartupMessages({
    library(CATALYST)
    library(cowplot)
    library(iCOBRA)
    library(ggplot2)
    library(methods)
    library(pkg)
})

# load data
data("sceEx", package="pkg")
data("idxEx", package="pkg")

# specify simulation parameters
n_cells <- c(800, 3e3)
n_genes <- 1250
n_samples <- 2
n_clusters <- 5
perc_type <- 0.2
perc_de <- c(0.01,0.3)
n_state <- (1-perc_type)*n_genes
n_tests <- n_clusters * n_state

methods <- pkg:::methods
diffcyt_methods <- c("limma", "LMM")
edgeR_methods <- c("rawCounts", "normCounts", "scaledCPM")

set.seed(1)
sd_range <- list(c(.5,1), c(.25,.5), c(.1,.25))
mean_fc <- c(2,1,.5)
mode_fc <- lapply(mean_fc, function(m) rnorm(1e2,m,.5))
for (i in seq_along(mode_fc)) 
    mode_fc[[i]] <- mode_fc[[i]][mode_fc[[i]] > 0]

# ------------------------------------------------------------------------------
cats <- c("DE", "DP", "DM", "DB")
for (cat in cats) {
    cat_probs <- as.integer(cats == cat)
    
    # initialise results
    pvals <- matrix(NA, nrow=n_tests, ncol=length(methods), dimnames=list(NULL, methods))
    if (cat == "DE") {
        pvals <- replicate(length(sd_range), pvals, simplify=FALSE)
        is_de <- matrix(NA, nrow=n_tests, ncol=length(sd_range))
    } else {
        pvals <- replicate(length(mode_fc), pvals, simplify=FALSE)
        is_de <- matrix(NA, nrow=n_tests, ncol=length(mode_fc))
    }
    
    for (i in seq_along(pvals)) {
        cat(paste0(i, "/", length(pvals), ":\n"))
        
        # simulate data
        if (cat == "DE") {
            sim <- simDD(sceEx, idxEx, 
                n_cells, n_genes, n_samples, n_clusters, 
                TRUE, TRUE, perc_type, perc_de, cat_probs, 
                sd_range=sd_range[[i]], seed=i)
        } else {
            sim <- simDD(sceEx, idxEx, 
                n_cells, n_genes, n_samples, n_clusters, 
                TRUE, TRUE, perc_type, perc_de, cat_probs, 
                mode_fc=mode_fc[[i]], seed=i)
        }
        
        # differential analysis
        for (j in 1:2) {
            pvals[[i]][,j] <- runDiffcyt(sim, data="normcounts", method=diffcyt_methods[j])$p_adj
            cat(".")
        }
        for (j in 1:3) {
            pvals[[i]][,j+2] <- runEdgeR(sim, method=edgeR_methods[j])$p_val
            cat(".")
        }
        
        # store gene information
        gi <- S4Vectors::metadata(sim)$gene_info
        gi <- gi[gi$gene %in% state_markers(sim), ]
        o <- c(sapply(seq_along(state_markers(sim)), function(i) 
            seq(i, n_clusters*n_state, n_state)))
        is_de[, i] <- gi$is_de[o]
        cat("\n")
    }
    saveRDS(pvals, file.path(path, sprintf("simDD_%s_pvals.rds", cat)))
    saveRDS(is_de, file.path(path, sprintf("simDD_%s_is_de.rds", cat)))
}

# ------------------------------------------------------------------------------
# plot results
for (cat in cats) {
    pvals <- readRDS(file.path(path, sprintf("simDD_%s_pvals.rds", cat)))
    is_de <- readRDS(file.path(path, sprintf("simDD_%s_is_de.rds", cat)))
    
    if (cat == "DE") {
        nms <- sapply(sd_range, function(i) sprintf("[%s,%s]", i[1], i[2]))
        nms <- sprintf("sd_range=%s", nms)
    } else {
        nms <- sprintf("mean(mode_fc)=%s", mean_fc)
    }
    ps <- lapply(seq_along(pvals), function(i) {
        cobradata <- COBRAData(
            pval=data.frame(pvals[[i]], check.names=FALSE),
            truth=data.frame(is_de=as.integer(is_de[, i])))
        cobradata <- calculate_adjp(cobradata)
        perf <- calculate_performance(cobradata,
            binary_truth="is_de", aspects=c("fdrtpr", "fdrtprcurve"))
        p <- plotPerf(perf, NULL, pkg:::method_colors) +
            annotate("text", x=1, y=0, label=nms[i], 
                size=2, hjust=1, vjust=0)
        if (i != 1)
            p <- p + labs(y=NULL)
        if (cat != "DB")
            p <- p + labs(x=NULL)
        return(p)
    })
    assign(sprintf("ps%s", cat), ps)
}
lgd <- get_legend(psDE[[1]] + 
        guides(color=guide_legend(override.aes=list(size=1))) +
        theme(legend.direction="horizontal"))
ps <- c(psDE, psDP, psDM, psDB)
for (i in seq_along(ps))
    ps[[i]] <- ps[[i]] + theme(legend.position="none")
labs <- rep("", length(ps))
labs[seq(1,length(ps),3)] <- cats
ps <- plot_grid(plotlist=ps, align="hv", axis="l",
    ncol=3, labels=labs, label_size=8, vjust=1, hjust=0)
ps <- plot_grid(ps, lgd, ncol=1, rel_heights=c(1,.02))
ggsave(file.path(path, "tpr_vs_fdr.pdf"), 
    plot=ps, width=12, height=18, unit="cm")
 