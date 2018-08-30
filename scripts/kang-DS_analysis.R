# Kang et al., DS analysis
# ==============================================================================
# - differential analysis using diffcyt & edgeR methods
# - comparison with published results
# ------------------------------------------------------------------------------

path <- "results/kang"

# load packages
suppressPackageStartupMessages({
    library(CATALYST)
    library(ComplexHeatmap)
    library(data.table)
    library(dplyr)
    library(ggplot2)
    library(pkg)
    library(RColorBrewer)
    library(reshape2)
    library(scater)
    library(SingleCellExperiment)  
})

# load data
daf <- readRDS(file.path(path, "kang_daf.rds"))
colData(daf)$marker_class <- "state" # make sure to check all genes

# initialise results
methods <- pkg:::methods
methods <- c(paste(rep(methods[1:2], each=2), c("meanNormcounts", "meanExprs"), sep="_"), methods[-c(1:2)])
n_clusters <- nlevels(cluster_ids(daf))
pvals <- matrix(NA,
    nrow=length(state_markers(daf))*n_clusters, ncol=length(methods),
    dimnames=list(NULL, methods))

# DE analysis
pvals[,1] <- runDiffcyt(daf, method="limma", data="normcounts", fun="mean")$p_adj
pvals[,2] <- runDiffcyt(daf, method="limma", data="exprs",      fun="mean")$p_adj

pvals[,3] <- runDiffcyt(daf, method="LMM", data="normcounts", fun="mean")$p_adj
pvals[,4] <- runDiffcyt(daf, method="LMM", data="exprs",      fun="mean")$p_adj

pvals[,5] <- runEdgeR(daf, method="rawCounts" )$p_val
pvals[,6] <- runEdgeR(daf, method="normCounts")$p_val
pvals[,7] <- runEdgeR(daf, method="scaledCPM" )$p_val

write.csv(pvals, file.path(path, "kang_pvals.csv"))

pvals <- read.csv(file.path(path, "kang_pvals.csv"), 
    header=TRUE, row.names=1, check.names=FALSE)

# reorder & split by cluster
pvals <- data.frame(pvals, check.names=FALSE)
pvals$cluster_id <- rep(levels(cluster_ids(daf)), length(state_markers(daf)))
o <- sapply(seq_len(n_clusters), function(i) seq(i, nrow(pvals), n_clusters))
pvals <- pvals[o, ]
rownames(pvals) <- NULL
pvals <- split(data.table(pvals), by="cluster_id", keep.by=FALSE)
pvals <- lapply(pvals, data.frame, check.names=FALSE)
for (i in seq_along(pvals))
    rownames(pvals[[i]]) <- state_markers(daf)

# get significant genes @ FDR 5%
s <- lapply(pvals, apply, 2, "<", 0.05)
s <- lapply(s, apply, 2, function(idx) {
    idx[is.na(idx)] <- FALSE
    state_markers(daf)[idx]
})

# compute expr. means by cluster
es <- assays(daf)$normcounts
es <- data.table(es, cluster_id=cluster_ids(daf))
es <- split(es, by="cluster_id", keep.by=FALSE)
means <- sapply(es, colMeans)
colnames(means) <- unique(cluster_ids(daf))
means <- means[, levels(cluster_ids(daf))]

# filter genes below mean 1
gs_keep <- apply(means, 2, function(gs) rownames(means)[gs >= 1])
s <- setNames(lapply(seq_along(s), function(k)
    setNames(lapply(seq_along(methods), function(m) {
        keep <- s[[k]][[m]] %in% gs_keep[[k]]
        s[[k]][[m]][keep]
    }), methods)), names(s))

# comparison w/ publication results
# ------------------------------------------------------------------------------

method_colors <- pkg:::method_colors
method_colors <- c(method_colors[1:2], c("purple3", "plum"), method_colors[-c(1:2)])
names(method_colors)[1:4] <- methods[1:4]
cluster_colors <- CATALYST:::cluster_cols[seq(1,20,2)[-c(8,9)]]

# read in cell-type specific DE genes
de_gs <- lapply(seq_len(n_clusters), function(i) {
    tbl <- readxl::read_xlsx("kang/de_genes.xlsx", sheet=i)
    tbl <- data.frame(tbl, row.names=1)
})
# reorder
names(de_gs) <- levels(cluster_ids(daf))[c(8,6,7,5,1,4,3,2)]
de_gs <- de_gs[levels(cluster_ids(daf))]

# plot nb. of DE genes by cluster
# & fraction of overlapping genes
df <- data.frame(t(sapply(s, sapply, length)), row.names=NULL, check.names=FALSE)
df <- df/sapply(de_gs, nrow)
df$cluster_id <- levels(cluster_ids(daf))
df <- melt(df, id.var="cluster_id")
df$value[df$value > 2] <- 2
p1 <- ggplot(df, aes(x=variable, y=value, fill=variable)) +
    geom_hline(yintercept=1, size=.2, lty=2) +
    geom_histogram(width=.8, stat="identity") +
    facet_wrap(~cluster_id, ncol=8, scale="free_x") +
    scale_x_discrete(NULL, expand=c(.2,0)) +
    scale_y_continuous("# DE genes (FDR = 5%)\nrelative to Kang et al.",
        limits=c(0,2), breaks=seq(0,2,.5), expand=c(.1,0)) +
    scale_fill_manual(NULL, values=method_colors) +
    theme_classic() + pkg:::prettify + theme(
        plot.margin=margin(0,.1,0,.1,"cm"),
        aspect.ratio=1,
        axis.title=element_text(size=5),
        strip.text=element_text(size=4, hjust=.5),
        panel.grid.major.y=element_line(size=.2, color="lightgrey"),
        axis.ticks.x=element_blank(), axis.text.x=element_blank(),
        legend.position="none", legend.direction="horizontal",
        legend.key.size=unit(.25,"cm"))

overlap <- sapply(seq_len(n_clusters), function(k)
    sapply(seq_along(methods), function(m) {
    gs <- s[[k]][[m]]
    x <- intersect(gs, rownames(de_gs[[k]]))
    length(x) / length(gs)
}))
rownames(overlap) <- methods
colnames(overlap) <- levels(cluster_ids(daf))
overlap <- melt(overlap)
overlap$value <- overlap$value
p2 <- ggplot(overlap, aes(x=Var1, y=value, fill=Var1)) +
    geom_histogram(width=.8, stat="identity") +
    facet_wrap(~Var2, ncol=8, scales="free_x") +
    scale_x_discrete(NULL, expand=c(.2,0)) +
    scale_y_continuous("Fraction of DE genes\nmatched with Kang et al.",
        limits=c(0,1), breaks=seq(0,1,.2), expand=c(.1,0)) +
    scale_fill_manual(NULL, values=method_colors) +
    theme_classic() + pkg:::prettify + theme(
        plot.margin=margin(0,.1,0,.1,"cm"),
        aspect.ratio=1,
        axis.title=element_text(size=5),
        strip.text=element_blank(),
        panel.grid.major.y=element_line(size=.2, color="lightgrey"),
        axis.ticks.x=element_blank(), axis.text.x=element_blank(),
        legend.position="bottom", legend.direction="horizontal",
        legend.key.size=unit(.25,"cm"))

ggsave(plot=p1, width=14, height=2.4, unit="cm",
    file.path(path, "kang_nb_de_gs.pdf"))
ggsave(plot=p2, width=14, height=3, unit="cm",
    file.path(path, "kang_overlap.pdf"))

# ------------------------------------------------------------------------------
# plot p-value distributions
df <- lapply(seq_along(pvals), function(i) pvals[[i]][rownames(de_gs[[i]]), ])
names(df) <- names(de_gs)
df <- lapply(df, melt)
df <- do.call(rbind, df)
df$cluster_id <- gsub(".[0-9]*$", "", rownames(df))
df <- df[!is.na(df$value), ]
df$variable <- gsub("_", "\n", df$variable)
df$variable <- gsub("meanNormcounts", "(meanNormcounts)", df$variable)
df$variable <- gsub("meanExprs", "(meanExprs)", df$variable)
df$variable <- factor(df$variable, levels=unique(df$variable))
p <- ggplot(df, aes(x=value, y=..scaled.., col=cluster_id)) +
    geom_density(size=.4) + 
    guides(color=guide_legend(override.aes=list(size=1))) +
    scale_color_manual(values=cluster_colors) +
    facet_wrap(~variable, ncol=4, scale="free") +
    labs(x="p-value", y="scaled density") +
    scale_x_continuous(limits=c(0,.6), breaks=seq(0,1,.2), expand=c(.1,0)) +
    scale_y_continuous(limits=c(0,1), breaks=seq(0,1,.2), expand=c(.15,0)) +
    theme_classic() + pkg:::prettify + theme(
        aspect.ratio=2/3, legend.key.size=unit(.25, "cm"),
        legend.position="bottom", legend.direction="horizontal",
        strip.text=element_text(size=6, hjust=.5),
        panel.grid.major=element_line(size=.2, color="lightgrey"))
ggsave(plot=p, width=14, height=7.6, unit="cm",
    file.path(path, "kang_pvals.pdf"))

# ------------------------------------------------------------------------------
# expr. heatmap for undetected genes w/ lowest p-values

n <- 10
diff <- sapply(seq_len(n_clusters), function(k) {
        inds <- !rownames(de_gs[[k]]) %in% unlist(s[[k]])
        top <- order(de_gs[[k]]$padj[inds])[seq_len(n)]
        rownames(de_gs[[k]])[top]
})
colnames(diff) <- levels(cluster_ids(daf))

es <- assays(daf)$exprs[, c(diff)]
df <- data.frame(es, cluster_id=cluster_ids(daf), sample_id=sample_ids(daf))
df <- df %>% group_by_("cluster_id", "sample_id") %>% summarise_all(mean)
df <- split(data.table(df), by="cluster_id", keep.by=FALSE)
for (i in seq_along(df)) {
    df[[i]] <- data.frame(df[[i]], row.names=1)
    df[[i]] <- t(df[[i]][, diff[, names(df)[i]]])
}
df <- df[levels(cluster_ids(daf))]
df <- do.call(rbind, df)

row_anno <- rowAnnotation(
    df=data.frame(cluster_id=factor(rep(levels(cluster_ids(daf)), each=n)), check.names=FALSE),
    col=list("cluster_id"=setNames(cluster_colors, levels(cluster_ids(daf)))),
    annotation_legend_param=list(
        title_gp=gpar(fontsize=8, fontface="bold", lineheight=.8),
        labels_gp=gpar(fontsize=6)))
hm <- Heatmap(
    mat=df, name="mean expression",
    col=rev(brewer.pal(11, "RdYlBu")),
    rect_gp=gpar(col="white"),
    column_order=levels(sample_ids(daf)),
    cluster_rows=FALSE, cluster_columns=FALSE,
    column_title="sample_id", column_title_side="bottom",
    split=rep(levels(cluster_ids(daf)), each=n), gap=unit(.5, "cm"),
    combined_name_fun=NULL,
    column_title_gp=gpar(fontsize=8, fontface="bold"),
    row_names_gp=gpar(fontsize=6),
    row_title_gp=gpar(fontsize=6),
    column_names_gp=gpar(fontsize=8),
    heatmap_legend_param=list(
        title_gp=gpar(fontsize=8, fontface="bold", lineheight=0.8),
        labels_gp=gpar(fontsize=6)))

pdf(file.path(path, "kang_top_undetected.pdf"), width=14/2.53, height=20/2.53)
row_anno + hm
decorate_heatmap_body("mean expression", {
    grid.lines(8/16, c(-8.95,.98), gp=gpar(lwd=1.5, lty=2, col="black"))
})
dev.off()

# ------------------------------------------------------------------------------
# expr. heatmap for Kang et al. DE genes w/ highest p-values
n <- 10
diff <- setNames(lapply(seq_len(n_clusters), function(k) { 
    setNames(lapply(seq_along(methods), function(m) {
        ps <- pvals[[k]][rownames(de_gs[[k]]), m]
        top <- order(ps, decreasing=TRUE, na.last=TRUE)[seq_len(n)]
        rownames(de_gs[[k]])[top]
    }), methods)
}), levels(cluster_ids(daf)))

exprs <- assays(daf)$exprs[, unique(unlist(diff))]
dfs <- setNames(lapply(levels(cluster_ids(daf)), function(k) {
    inds <- cluster_ids(daf) == k
    es <- exprs[inds, unlist(diff[[k]])]
    df <- data.frame(es, sample_id=sample_ids(daf)[inds])
    df <- df %>% group_by_("sample_id") %>% summarise_all(mean)
    df <- t(data.frame(df, row.names=1))
    rownames(df) <- gsub("\\.+[0-9]", "", rownames(df))
    return(df)
}), levels(cluster_ids(daf)))

split <- anno_df <- factor(rep(methods, each=n), levels=methods)
for (k in levels(cluster_ids(daf))) {
    pdf(file.path(path, sprintf("kang_highest_pvals_%s.pdf", k)),
        width=14/2.53, height=20/2.53)
    row_anno <- rowAnnotation(
        df=data.frame(method=anno_df, check.names=FALSE),
        col=list("method"=method_colors),
        annotation_legend_param=list(
            title_gp=gpar(fontsize=8, fontface="bold", lineheight=.8),
            labels_gp=gpar(fontsize=6)))
    hm <- Heatmap(
        mat=dfs[[k]], name="mean expression",
        col=rev(brewer.pal(11, "RdYlBu")),
        rect_gp=gpar(col="white"),
        column_order=levels(sample_ids(daf)),
        cluster_rows=FALSE, 
        cluster_columns=FALSE,
        column_title=k,
        split=split, gap=unit(.5, "cm"),
        column_title_gp=gpar(fontsize=8, fontface="bold"),
        combined_name_fun=NULL,
        row_title_gp=gpar(fontsize=6),
        row_names_gp=gpar(fontsize=6),
        column_names_gp=gpar(fontsize=8),
        heatmap_legend_param=list(
            title_gp=gpar(fontsize=8, fontface="bold", lineheight=.8),
            labels_gp=gpar(fontsize=6)))
    print(row_anno + hm)
    decorate_heatmap_body("mean expression", {
        grid.lines(8/16, c(-7.43,.99), gp=gpar(lwd=1.5, lty=2, col="black"))
    })
    dev.off()
}

