# simDD, Seurat clustering
# ==============================================================================
# This script investigates how well clusters 
# simulated by simDD can be retrieved by Seurat.
# ------------------------------------------------------------------------------

path <- "results/simDD-seurat"

# load packages
suppressPackageStartupMessages({
    library(CATALYST)
    library(dplyr)
    library(ggplot2)
    library(methods)
    library(pkg)
    library(reshape2)
    library(Seurat)
    library(SummarizedExperiment)
    library(zeallot)
})

# load data
data("sceEx", package="pkg")
data("idxEx", package="pkg")

# specify simulation parameters
n_sims <- 10

# initialise results
cluster_ids <- type_markers <- gis <- vector("list", n_sims)

for (i in seq_len(n_sims)) {
    cat(i, "/", n_sims, "\n", sep="")
    
    set.seed(i)
    n_clusters <- sample(4:12,1)
    sd_range <- sort(sample(runif(2,1,3),2))
    mode_fc <- rnorm(1e2,runif(1,.5,4),.5)
    mode_fc <- mode_fc[mode_fc > 0]
    cat_probs <- sample(runif(4))
    cat_probs <- cat_probs/sum(cat_probs)
    
    # simulate data
    sim <- simDD(sceEx, idxEx, n_cells=c(100,1e3), n_genes=1250,
        n_clusters=n_clusters, n_samples=2, rand_clusters=TRUE, rand_samples=TRUE,
        sd_range=sd_range, mode_fc=mode_fc,
        perc_type=.2, perc_de=c(0,1), cat_probs=cat_probs, seed=i)
    
    # Seurat clustering
    message("Clustering...")
    counts <- t(assays(sim)$counts)
    colnames(counts) <- sprintf("cell%s", seq_len(ncol(counts)))
    so <- CreateSeuratObject(counts)
    so <- NormalizeData(so, display.progress=FALSE)
    so <- FindVariableGenes(so, display.progress=FALSE, do.plot=FALSE)
    so <- ScaleData(so, display.progress=FALSE)
    so <- RunPCA(so, do.print=FALSE)
    so <- FindClusters(so, so@var.genes, "pca", 1:10, print.output=FALSE)
    
    # store gene information
    n_clusters <- nlevels(cluster_ids(sim))
    n_state <- length(state_markers(sim))
    gi <- metadata(sim)$gene_info
    gi <- gi[gi$gene %in% state_markers(sim), ]
    o <- c(sapply(seq_along(state_markers(sim)), function(i) 
        seq(i, n_clusters*n_state, n_state)))
    gis[[i]] <- gi[o, ]
    
    # get genes used for clustering & cluster IDs
    cluster_ids[[i]] <- data.frame(
        truth=cluster_ids(sim), 
        seurat=factor(as.numeric(as.character(so@ident)) + 1))
    type_markers[[i]] <- list(
        truth=type_markers(sim),
        seurat=so@var.genes)
    cat("\n")
}

saveRDS(cluster_ids, file.path(path, "seurat_cluster_ids.rds"))
saveRDS(type_markers, file.path(path, "seurat_cluster_gs.rds"))
saveRDS(gis, file.path(path, "seurat_gis.rds"))

# get confusion matrix (nb. of TP, FP, TN, FN)
source("utils.R")
cm <- lapply(cluster_ids, function(ids) 
    get_status(ids[, 1], ids[, 2]))

# compute evaluation metrics
funs <- c("precision", "accuracy", "recall", "f1")
res <- sapply(funs, function(fun) sapply(cm, fun))
res <- melt(data.frame(res))

# rename & reorder
res$variable <- as.character(res$variable)
res$variable[res$variable == "f1"] <- "F1 score"
res$variable <- factor(res$variable, levels=c(
    "accuracy", "precision", "recall", "F1 score"))

p <- ggplot(res, aes(y=value, col=variable)) + 
    guides(color=guide_legend(NULL)) +
    geom_boxplot(width=.8, size=.4) +
    scale_x_discrete(NULL, expand=c(.15,0)) +
    scale_y_continuous(NULL, limits=c(.2,1), breaks=seq(.2,1,.2), expand=c(.1,0)) +
    theme_classic() + pkg:::prettify + theme(
        aspect.ratio=2/3,
        panel.grid.major=element_line(size=.2, color="lightgrey"),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        legend.key.size=unit(.25, "cm"))
ggsave(file.path(path, "simDD-seurat_scores.pdf"),
    plot=p, width=7, height=3.2, unit="cm")

