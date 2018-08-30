# simDD, null simulation
# ==============================================================================
# Generates p-value distributions for 3 null stimulation replicates
# (no differentially expressed genes, 1000 tests per run)
# ------------------------------------------------------------------------------

path <- "results/simDD-null_sim"

# load packages
suppressPackageStartupMessages({
    library(CATALYST)
    library(pkg)
    library(ggplot2)
    library(methods)
})

# load data
data("sceEx", package="pkg")
data("idxEx", package="pkg")

# simulate data
n_reps <- 3
sim <- lapply(seq_len(n_reps), function(i)
    simDD(sceEx, idxEx, 
        n_cells=c(2e3, 4e3), n_genes=250, n_clusters=4, n_samples=4, 
        rand_clusters=TRUE, rand_samples=TRUE, perc_type=0, perc_de=0, seed=i))

# initialise results
methods <- pkg:::methods
methods <- c(methods[1:2], "FDA", methods[-c(1:2)])
n_gs <- length(state_markers(sim[[1]]))
n_cs <- nlevels(cluster_ids(sim[[1]]))
pvals <- replicate(n_reps, simplify=FALSE, 
    matrix(NA, nrow=n_gs * n_cs, ncol=length(methods), 
        dimnames=list(rep(state_markers(sim[[1]]), n_cs), methods)))

# differential analysis
for (i in seq_len(n_reps)) {
    pvals[[i]][, 1] <- runDiffcyt(sim[[i]], method="limma", data="normcounts", log=FALSE, fun="mean")$p_val; cat(".")
    pvals[[i]][, 2] <- runDiffcyt(sim[[i]], method="LMM"  , data="normcounts", log=FALSE, fun="mean")$p_val; cat(".")
    
    pvals[[i]][, 3] <- runFDA(sim[[i]], reso=25, n_perm=50)$p_val; cat(".")
    
    pvals[[i]][, 4] <- runEdgeR(sim[[i]], method="rawCounts" )$p_val; cat(".")
    pvals[[i]][, 5] <- runEdgeR(sim[[i]], method="normCounts")$p_val; cat(".")
    pvals[[i]][, 6] <- runEdgeR(sim[[i]], method="scaledCPM" )$p_val; cat(".")
    
    cat("\n")
}
saveRDS(pvals, file.path(path, "pvals.rds"))

# plot results
df <- data.frame(do.call(rbind, pvals), check.names=FALSE, row.names=NULL)
df$seed <- factor(rep(seq_len(n_reps), each=nrow(pvals[[1]])))
df <- reshape2::melt(df, id.var="seed")

p <- ggplot(df, aes(x=value, y=..scaled.., lty=seed)) + 
    facet_wrap(~ variable, scale="free", nrow=2) +
    geom_density(alpha=0.2, adjust=0.5, col="darkblue", fill="royalblue") +
    scale_linetype_manual(values=1:3) +
    scale_x_continuous(limits=c(0,1), breaks=seq(0,1,.25), expand=c(.05,0)) +
    scale_y_continuous(limits=c(0,1), breaks=seq(0,1,.2), expand=c(.075,0)) +
    labs(x="p-value", y="scaled density") +
    theme_classic() + pkg:::prettify + theme(
        aspect.ratio=2/3,
        legend.margin=margin(0,0,0,0, "cm"),
        strip.text=element_text(hjust=0.5))
ggsave(file.path(path, "simDD-null_sim.pdf"), 
    plot=p, width=14, height=7.3, unit="cm")

