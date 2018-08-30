# Exemplary ECDFs
# ==============================================================================
# Plots an exemplary set of ECDFs as computed by FDA 
# for 2 genes per categories.
# ------------------------------------------------------------------------------

# load packages
suppressPackageStartupMessages({
    library(CATALYST)
    library(data.table)
    library(ggplot2)
    library(pkg)
    library(reshape2)
    library(SummarizedExperiment)
})

# load data
data("sceEx", package="pkg")
data("idxEx", package="pkg")

# simulate data
sim <- simDD(sceEx, idxEx, 
    n_cells=6e3, n_genes=500, n_clusters=1, n_samples=3, 
    perc_type=0, perc_de=.5, cat_probs=rep(.25, 4), seed=42)

# subset 2 genes per category for plotting
gi <- metadata(sim)$gene_info
set.seed(3)
genes_keep <- c(sapply(levels(gi$cat), function(c) sample(which(gi$cat == c), 2)))

# split exprs. by sample
es <- assays(sim)$exprs[, genes_keep]
dt <- data.table(es, sample_id=sample_ids(sim))
es_by_sample <- split(dt, by="sample_id", keep.by=FALSE)

# get grid-points to evaluate ECDFs at
xs <- apply(colRanges(es), 1, function(rng) seq(rng[1], rng[2], length=50))

# get ECDFs by sample & gene
ecdfs <- lapply(es_by_sample, sapply, stats::ecdf)
ecdfs.x <- lapply(ecdfs, function(ecdf) 
    sapply(seq_along(genes_keep), function(g) ecdf[[g]](xs[, g])))

# plot ECDFs
df <- do.call(rbind, ecdfs.x)
colnames(df) <- paste0(gi$gene, sprintf("(%s)", gi$cat))[genes_keep]
df <- data.frame(check.names=FALSE, xs, df,
    sample_id=rep(names(ecdfs), each=nrow(xs)))
colnames(df)[seq_len(ncol(xs))] <- "x"
df$condition <- gsub("[0-9]", "", df$sample_id)
df <- melt(df, id.var=c("condition", "sample_id", "x"))

# reorder
o <- c(9,10,11,12,1,4,6,7,2,3,5,8)
df$variable <- factor(df$variable, levels=levels(df$variable)[o])

p <- ggplot(df, aes(x=x, y=value, group=sample_id, col=condition)) +
    facet_wrap(~variable, ncol=4, scales="free") + geom_path(size=.4) + 
    scale_x_continuous(limits=c(0,10), breaks=seq(0,10,2.5), expand=c(.1,0)) +
    scale_y_continuous(limits=c(0,1), breaks=seq(0,1,.25), expand=c(.15,0)) +
    labs(x="expression", y="cumulative probability") + 
    theme_classic() + pkg:::prettify + theme(
        aspect.ratio=2/3, legend.position="none",
        strip.text=element_text(hjust=0.5),
        panel.grid.major=element_line(size=.2, color="lightgrey"))

ggsave("figures/ECDFs.pdf", 
    plot=p, width=14, height=9.2, unit="cm")
