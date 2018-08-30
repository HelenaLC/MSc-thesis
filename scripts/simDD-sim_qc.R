# simDD, quality control
# ==============================================================================
# Generates basic plots for quality control,
# comparing simDD simulated data & data from Koh et al.
# - variance explained by cluster ID, sample ID & condition factors
# - expression frequency vs. mean
# - top expressed genes
# - distribution of library sizes
# - dispersion vs. mean
# ------------------------------------------------------------------------------

# load packages
suppressPackageStartupMessages({
    library(cowplot)
    library(DESeq2)
    library(edgeR)
    library(pkg)
    library(MultiAssayExperiment)
    library(scater)
    library(SingleCellExperiment)
    library(SummarizedExperiment)
})

# load data
data("idxEx", package="pkg")
data("sceEx", package="pkg")

# simulate data
sim <- simDD(sceEx, idxEx, 
    n_cells=800, n_genes=1250, n_samples=2, n_clusters=3,
    perc_type=.2, perc_de=.05, cat_probs=c(.7,.2,.08,.02), seed=22)

# load reference 
koh <- readRDS("koh/SRP073808.rds")
cts <- assays(experiments(koh)[["gene"]])[["count"]]
tpm <- assays(experiments(koh)[["gene"]])[["TPM"]]

# prep. data
sim <- SingleCellExperiment(
    assays=lapply(assays(sim),t),
    colData=rowData(sim))
koh <- SingleCellExperiment(
    assays=list(counts=cts, tpm=tpm),
    colData=colData(koh))
koh <- koh[rowSums(assay(koh) > 0) != 0, ]
rownames(koh) <- gsub("\\.+[0-9]*", "", rownames(koh))

# calculate exprs.
assays(sim)$logcounts <- log2(assays(sim)$normcounts+1)
assays(koh)$logcounts <- log2(calculateCPM(koh)+1)

# calculate QC metrics
suppressMessages({
    koh <- calculateQCMetrics(koh)
    sim <- calculateQCMetrics(sim)
})

# ------------------------------------------------------------------------------
# variance explained
vars <- c("cluster_id", "condition", "sample_id")
p <- plotExplanatoryVariables(sim, variables=vars) +
    scale_y_continuous(limits=c(0,.85), expand=c(.075,0)) +
    labs(y="density") + theme_bw() + pkg:::prettify + 
    theme(aspect.ratio=2/3, legend.key.size=unit(.25,"cm"),
        panel.grid.major=element_line(size=.2, color="lightgrey"))
ggsave("results/simDD-sim_qc/qc_var_explained.pdf", 
    plot=p, width=8, height=4.2, unit="cm")

# ------------------------------------------------------------------------------
# expr. frequency vs. mean
source("scripts/plotExprsFreqVsMean2.R")
p1 <- plotExprsFreqVsMean2(sim) + pkg:::prettify + ggtitle("simDD") +
    scale_x_continuous(limits=c(0,15), breaks=seq(0,15,5), expand=c(.1,0)) +
    scale_y_continuous("expressing cells (%)", 
        limits=c(0,100), breaks=seq(0,100,20), expand=c(.1,0))
p2 <- plotExprsFreqVsMean2(koh) + labs(title="Koh et al.", y=NULL) + pkg:::prettify +
    scale_x_continuous(limits=c(0,15), breaks=seq(0,15,5), expand=c(.1,0)) +
    scale_y_continuous(limits=c(0,100), breaks=seq(0,100,20), expand=c(.1,0))
p <- plot_grid(p1, p2, align="v")
ggsave("results/simDD-sim_qc/qc_expr_freq_vs_mean.pdf",
    plot=p, width=12, height=6, unit="cm")

# ------------------------------------------------------------------------------
# top expressed genes
p1 <- plotHighestExprs(sim, n=25, colour_cells_by=NULL) +
    guides(fill=FALSE) + labs(title="simDD", y=NULL) + 
    theme_bw() + pkg:::prettify +
    theme(axis.text.y=element_text(size=5), plot.title=element_text(size=7))
p2 <- plotHighestExprs(koh, n=25, colour_cells_by=NULL) +
    guides(fill=FALSE) + labs(title="Koh et al.", y=NULL) + 
    theme_bw() + pkg:::prettify +
    theme(axis.text.y=element_text(size=5), plot.title=element_text(size=7))
p <- plot_grid(p1, p2, align="v")
ggsave(plot=p, "results/simDD-sim_qc/qc_top_expr.pdf", 
    width=12, height=5.6, unit="cm")

# ------------------------------------------------------------------------------
# distribution of library sizes
lib_sizes_sim <- colSums(assay(sim))
lib_sizes_koh <- colSums(assay(koh))
df_sim <- data.frame(x=lib_sizes_sim)
df_koh <- data.frame(x=lib_sizes_koh)
p1 <- ggplot(df_sim, aes(x=x, y=..ncount..)) + 
    labs(title="simDD", x="total library size", y="normalised count") +
    scale_x_continuous(expand=c(.025,0),
        labels=function(x) format(x, scientific=TRUE)) +
    scale_y_continuous(limits=c(0,1), breaks=seq(0,1,.2), expand=c(.1,0)) + 
    geom_histogram(color="white", fill="blue3", bins=50) +
    theme_classic() + pkg:::prettify + theme(aspect.ratio=1/4,
        axis.text.x=element_text(angle=30, hjust=1, vjust=1))
p2 <- ggplot(df_koh, aes(x=x, y=..ncount..)) + 
    labs(title="Koh et al.", x="total library size", y=NULL) +
    scale_x_continuous(expand=c(.025,0)) +
    scale_y_continuous(limits=c(0,1), breaks=seq(0,1,.2), expand=c(.1,0)) +
    geom_histogram(color="white", fill="blue3", bins=50) +
    theme_classic() + pkg:::prettify + theme(aspect.ratio=1/4,
        axis.text.x=element_text(angle=30, hjust=1, vjust=1))
p <- plot_grid(p1, p2, align="v")
ggsave("results/simDD-sim_qc/qc_lib_sizes.pdf", 
    plot=p, width=14, height=3.4, unit="cm")

# ------------------------------------------------------------------------------
# dispersion vs. mean of normalised counts / TPM
means_sim <- rowMeans(assays(sim)$normcounts)
means_koh <- rowMeans(assays(koh)$tpm)
means_koh <- means_koh[means_koh > 1e-3]
disps_sim <- edgeR::estimateDisp(assay(sim), model.matrix(~colData(sim)$condition))
disps_koh <- edgeR::estimateDisp(assay(koh)[names(means_koh), ])
df_sim <- data.frame(x=means_sim, y=disps_sim$tagwise.dispersion, z=disps_sim$trended.dispersion)
df_koh <- data.frame(x=means_koh, y=disps_koh$tagwise.dispersion, z=disps_koh$trended.dispersion)
p1 <- ggplot() + 
    geom_point(data=df_sim, aes(x=x, y=y), size=.1, alpha=.4) + 
    geom_smooth(data=df_sim, aes(x=x, y=z), formula=y~z, size=.8, color="red2") +
    labs(title="simDD", x="mean normalised count", y="dispersion") +
    scale_x_continuous(trans="log10", limits=c(10, 1e4), expand=c(.1,0), 
        labels=function(x) format(x, scientific=TRUE)) +
    scale_y_continuous(trans="log10", limits=c(1e-1, 1e2), expand=c(.1,0)) + 
    theme_bw() + pkg:::prettify + theme(
        panel.grid.major=element_line(size=.2, color="lightgrey"))
p2 <- ggplot() + 
    geom_point(data=df_koh, aes(x=x, y=y), size=.1, alpha=.4) + 
    geom_smooth(data=df_koh, aes(x=x, y=z), formula=y~z, size=.8, color="red2") +
    labs(title="Koh et al.", x="mean TPM", y=NULL) +
    scale_x_continuous(trans="log10", limits=c(1e-3, 1e4), expand=c(.1,0)) +
    scale_y_continuous(trans="log10", limits=c(1e-3, 1e2), expand=c(.1,0)) + 
    theme_bw() + pkg:::prettify + theme(
        panel.grid.major=element_line(size=.2, color="lightgrey"))
p <- plot_grid(p1, p2, align="v")
ggsave("results/simDD-sim_qc/qc_disp_vs_mean.pdf", 
    plot=p, width=12, height=6.2, unit="cm")
