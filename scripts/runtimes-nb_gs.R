# runtimes for increasing nb. of genes
# ==============================================================================
# This script investigates the effect the nb. of tested genes has on the 
# runtimes of all differential analysis methods.
# ------------------------------------------------------------------------------

path <- "results/runtimes"

# load packages
suppressPackageStartupMessages({
    library(CATALYST)
    library(ggplot2)
    library(methods)
    library(pkg)
    library(reshape2)
})

# load data
data("sceEx", package="pkg")
data("idxEx", package="pkg")

# set simulation parameters
methods <- c(pkg:::methods, "FDA")
n_genes <- c(100, 250, 500, 1e3, 2e3, 5e3, 1e4)

# initialise results
rts <- matrix(NA,
    nrow=length(n_genes), 
    ncol=length(methods), 
    dimnames=list(n_genes, methods))

for (i in seq_along(n_genes)) {
    cat(n_genes[i], "\n")
    
    # simulate data
    sim <- simDD(sceEx, idxEx, 
        n_cells=2e3, n_genes[i], n_samples=4, n_clusters=1, 
        perc_type=0, perc_de=.2, cat_probs=c(1,0,0,0), seed=i)
    
    # DE analysis
    rts[i,1] <- system.time(runDiffcyt(sim, method="limma", data="normcounts"))[[3]]; cat(".")
    rts[i,2] <- system.time(runDiffcyt(sim, method="LMM"  , data="normcounts"))[[3]]; cat(".")
    
    rts[i,3] <- system.time(runEdgeR(sim, method="rawCounts") )[[3]]; cat(".")
    rts[i,4] <- system.time(runEdgeR(sim, method="normCounts"))[[3]]; cat(".")
    rts[i,5] <- system.time(runEdgeR(sim, method="scaledCPM") )[[3]]; cat(".")
    
    if (i < 6)
        rts[i,6] <- system.time(runFDA(sim, n_perm=50))[[3]]; cat(".")
    cat("\n")
}  
saveRDS(rts,   file.path(path, "runtimes-nb_gs.rds"))

# plot runtimes vs. nb. of genes
rts <- readRDS(file.path(path, "runtimes-nb_gs.rds"))
rts <- melt(rts)
rts <- rts[!is.na(rts$value), ]
p <- ggplot(rts, aes(x=Var1, y=value, col=Var2)) + 
    geom_point(size=1.5) + geom_path(lty=2, size=.5) +
    scale_x_continuous("n_genes", trans="log10", 
        limits=c(100,1e4), breaks=n_genes[-c(1,3)], expand=c(.1,0)) +
    scale_y_continuous("elapsed time (s)", trans="log10", 
        limits=c(1e-1,1e3), expand=c(.1,0)) +
    scale_color_manual(NULL, values=c(pkg:::method_colors, FDA="limegreen")) +
    guides(color=guide_legend(override.aes=list(size=1, lty=NA))) +
    theme_classic() + pkg:::prettify + theme(aspect.ratio=2/3,
        legend.key.size=unit(.25, "cm"),
        axis.text.x=element_text(angle=30, hjust=1, vjust=1),
        panel.grid.major.y=element_line(size=.2, color="lightgrey"))
ggsave(file.path(path, "runtimes-nb_gs.pdf"), 
    plot=p, width=10, height=5.2, unit="cm")

