# runtimes FDA n_perm
# ==============================================================================
# This script investigates the runtime of differential analysis with FDA
# for varying nb. of permutations evaluated by fda::tperm.fd.
# ------------------------------------------------------------------------------

path <- "results/runtimes"

# load packages
suppressPackageStartupMessages({
    library(ggplot2)
    library(matrixStats)
    library(methods)
    library(pkg)
    library(reshape2)
})

# load data
data("sceEx", package="pkg")
data("idxEx", package="pkg")

# specify simulation parameters
n_sims <- 10
n_perm <- seq(50, 500, 50)
rts <- matrix(NA, 
    nrow=n_sims, ncol=length(n_perm), 
    dimnames=list(NULL, n_perm))

for (i in seq_len(n_sims)) {
    cat(i, "/", n_sims, "\n", sep="")
    # simulate data
    sim <- simDD(sceEx, idxEx,  
        n_cells=3e3, n_genes=800, n_clusters=1, n_samples=3, 
        perc_type=.25, perc_de=.5, cat_probs=c(1,0,0,0), seed=i)
    for (j in seq_along(n_perm)) {
        # differential analysis
        rts[i, j] <- system.time(runFDA(sim, reso=50, n_perm=n_perm[j]))[[3]]
        cat(".")
    }  
    cat("\n")
}
write.csv(rts, file.path(path, "runtimes_FDA_nperm.csv"))

# plot results
p <- ggplot(melt(rts), aes(x=variable, y=value)) + 
    geom_boxplot(size=.2, width=.8, fill="aliceblue", color="blue3", outlier.color=NA) +
    scale_x_discrete(expand=c(.1,0)) +
    scale_y_continuous(limits=c(0,60), breaks=seq(0,60,20), expand=c(.15,0)) +
    labs(x="n_perm", y="elapsed time (s)") + theme_classic() + pkg:::prettify +
    theme(aspect.ratio=2/3,
        panel.grid.major=element_line(size=.2, color="lightgrey"))

ggsave(plot=p, width=7, height=5, unit="cm",
    file.path(path, "runtimes_FDA_nperm.pdf"))

