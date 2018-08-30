# runtimes FDA reso
# ==============================================================================
# This script investigates the runtime of differential analysis with FDA
# for varying nb. of grid-points evaluated by fda::tperm.fd.
# ------------------------------------------------------------------------------

path <- "results/runtimes"

# load packages
suppressPackageStartupMessages({
    library(ggplot2)
    library(methods)
    library(pkg)
    library(reshape2)  
})

# load data
data("sceEx", package="pkg")
data("idxEx", package="pkg")

# specify simulation parameters
n_sims <- 10
reso <- seq(10, 100, 10)
rts <- matrix(NA, 
    nrow=n_sims, ncol=length(reso), 
    dimnames=list(NULL, reso))

for (i in seq_len(n_sims)) {
    cat(i, "/", n_sims, "\n", sep="")
    # simulate data
    sim <- simDD(sceEx, idxEx,  
        n_cells=2e3, n_genes=500, n_clusters=1, n_samples=2, 
        perc_type=0, perc_de=.5, cat_probs=c(1,0,0,0), seed=i)
    for (j in seq_along(reso)) {
        # differential analysis
        rts[i, j] <- system.time(runFDA(sim, reso=reso[j], n_perm=200))[[3]]
        cat(".")
    }  
    cat("\n")
}
write.csv(rts, file.path(path, "runtimes_FDA_reso.csv"))

# plot results
p <- ggplot(melt(rts), aes(x=variable, y=value)) + 
    geom_boxplot(size=.2, width=.8, fill="aliceblue", color="blue3", outlier.color=NA) +
    scale_x_discrete(expand=c(.1,0)) +
    scale_y_continuous(limits=c(0,45), breaks=seq(0,45,15), expand=c(.15,0)) +
    labs(x="reso", y="elapsed time (s)") + theme_classic() + pkg:::prettify +
    theme(aspect.ratio=2/3,
        panel.grid.major=element_line(size=.2, color="lightgrey"))

ggsave(plot=p, width=7, height=5, unit="cm",
    file.path(path, "runtimes_FDA_reso.pdf"))
