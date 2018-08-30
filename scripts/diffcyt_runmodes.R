# diffcyt performance for 6 different runmodes
# ==============================================================================
# This script investigates the performance of diffcyt for varying
# - data inputs (raw counts, normalised counts, expressions) &
# - summary statistics (means, medians).
# ------------------------------------------------------------------------------

path <- "results/diffcyt_runmodes"

# load packages
suppressPackageStartupMessages({
    library(CATALYST)
    library(cowplot)
    library(ggplot2)
    library(iCOBRA)
    library(methods)
    library(pkg)
})

# load data
data("sceEx", package="pkg")
data("idxEx", package="pkg")

methods <- c("meanRawCounts", "medRawCounts", "meanNormCounts", "medNormCounts",  "meanExprs", "medExprs")
gene_cats <- c("DE", "DP", "DM", "DB")
cat_probs <- matrix(c(diag(4), rep(.25,4)), byrow=TRUE, nrow=5, 
    dimnames=list(c(gene_cats, "mixed"), gene_cats))

# specify simulation parameters
n_cells <- c(2e3, 4e3)
n_genes <- 1250
n_samples <- 2
n_clusters <- 3
perc_type <- 0.2
perc_de <- c(.05,.3)
n_state <- (1-perc_type)*n_genes
n_tests <- n_clusters * n_state

# initialise results
res <- setNames(vector("list", nrow(cat_probs)), rownames(cat_probs))
is_de <- matrix(NA, nrow=n_tests, ncol=nrow(cat_probs))

for (i in seq_len(nrow(cat_probs))) {
    cat(rep("=", 10), " ", rownames(cat_probs)[i], " ", rep("=", 10), "\n", sep="")
    
    # simulate data
    sim <- simDD(sceEx, idxEx, 
        n_cells, n_genes, n_samples, n_clusters, TRUE, TRUE,
        perc_type, perc_de, cat_probs[i, ], seed=i)
    
    # initialise results
    pvals <- matrix(NA, nrow=n_tests, ncol=length(methods), 
        dimnames=list(rep(state_markers(sim), each=n_clusters), methods))
    
    # differential analysis
    pvals[, 1] <- runDiffcyt(sim, method="limma", data="counts",     log=FALSE, fun="mean")  $p_adj; cat(".")
    pvals[, 2] <- runDiffcyt(sim, method="limma", data="counts",     log=FALSE, fun="median")$p_adj; cat(".")
    pvals[, 3] <- runDiffcyt(sim, method="limma", data="normcounts", log=FALSE, fun="mean")  $p_adj; cat(".")
    pvals[, 4] <- runDiffcyt(sim, method="limma", data="normcounts", log=FALSE, fun="median")$p_adj; cat(".")
    pvals[, 5] <- runDiffcyt(sim, method="limma", data="normcounts", log=TRUE,  fun="mean")  $p_adj; cat(".")
    pvals[, 6] <- runDiffcyt(sim, method="limma", data="normcounts", log=TRUE,  fun="median")$p_adj; cat(".")
    
    # store indices of differentially expressed genes
    gi <- S4Vectors::metadata(sim)$gene_info
    gi <- gi[gi$gene %in% state_markers(sim), ]
    o <- c(sapply(seq_along(state_markers(sim)), function(i) 
        seq(i, n_clusters*n_state, n_state)))
    is_de[, i] <- gi$is_de[o]
    res[[i]] <- pvals
    cat("\n")
}
saveRDS(res,   file.path(path, "diffcyt_runmodes_pvals.rds"))
saveRDS(is_de, file.path(path, "diffcyt_runmodes_is_de.rds"))

pvals <- readRDS(file.path(path, "diffcyt_runmodes_pvals.rds"))
is_de <- readRDS(file.path(path, "diffcyt_runmodes_is_de.rds"))

# plot TPR vs. FDR curves
cols <- c("darkblue", "royalblue", "red2", "orange", "purple2", "plum") 
cols <- CATALYST:::cluster_cols[c(1:4,11:12)]
cols <- setNames(cols, methods)
ps <- vector("list", length(pvals))
for (i in seq_along(pvals)) {
    cobradata <- COBRAData(
        pval=data.frame(pvals[[i]], check.names=FALSE, row.names=NULL),
        truth=data.frame(is_de=as.integer(is_de[, i])))
    cobradata <- calculate_adjp(cobradata)
    perf <- calculate_performance(cobradata, 
        binary_truth="is_de", aspects=c("fdrtpr", "fdrtprcurve"))
    ps[[i]] <- plotPerf(perf, rownames(cat_probs)[i], cols) +
        guides(color=guide_legend(override.aes=list(size=1))) +
        theme(legend.key.size=unit(.25, "cm"))
    if (i != 4 & i != 5) ps[[i]] <- ps[[i]] + xlab(NULL)
    if (i != 1 & i != 4) ps[[i]] <- ps[[i]] + ylab(NULL)
}

# arrange plots
lgd <- get_legend(ps[[1]])
for (i in seq_along(ps))
    ps[[i]] <- ps[[i]] + theme(legend.position="none")
ps <- align_plots(plotlist=ps, align="hv", axis="l")
top <- plot_grid(plotlist=ps[1:3], nrow=1)
bot <- plot_grid(ps[[4]], ps[[5]], lgd, nrow=1)
p <- plot_grid(top, bot, ncol=1)
ggsave(plot=p, width=14, height=10.2, unit="cm",
    file.path(path, "diffcyt_runmodes_tpr_vs_fdr.pdf"))

