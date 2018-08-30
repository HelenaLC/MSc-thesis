# Kang et al., data overview
# ==============================================================================
# Generates plots for general data overview
# - cell-type compositions by sample
# - t-SNE projection
# - # of cells by cluster & sample
# ------------------------------------------------------------------------------

path <- "results/kang"

# load packages
library(CATALYST)
library(cowplot)
library(ggplot2)

# load data
daf <- readRDS(file.path(path, "kang_daf.rds"))

# ------------------------------------------------------------------------------
# relative cluster abundances
p <- plotAbundances(daf) + pkg:::prettify + 
    theme(legend.key.size=unit(.25,"cm"))
ggsave(file.path(path, "kang_cluster_props.pdf"), 
    plot=p, width=12, height=5.6, unit="cm")

# ------------------------------------------------------------------------------
# t-SNE
df <- data.frame(
    rowData(daf)[metadata(daf)$tsne_inds, ],
    tSNE1=metadata(daf)$tsne$Y[, 1], 
    tSNE2=metadata(daf)$tsne$Y[, 2])

thm <- theme(legend.key.height=unit(.2,"cm"),
    panel.grid.major=element_line(size=.2, color="lightgrey"),
    axis.line=element_blank(), axis.title=element_text(face="bold"))

p1 <- ggplot(df, aes(x=tSNE1, y=tSNE2)) +
    geom_point(data=df, size=.75, aes(col=condition)) +
    guides(color=guide_legend(override.aes=list(size=2), ncol=1)) + 
    scale_color_manual(values=c("mediumblue", "limegreen")) +
    pkg:::prettify + thm

p2 <- ggplot(df, aes(x=tSNE1, y=tSNE2)) + ylab(NULL) +
    geom_point(data=df, size=.75, aes(col=cluster_id)) +
    guides(color=guide_legend(override.aes=list(size=2), ncol=1)) + 
    scale_color_manual(values=CATALYST:::cluster_cols[1:8]) +
    pkg:::prettify + thm

l1 <- get_legend(p1)
l2 <- get_legend(p2)
p1 <- p1 + theme(legend.position="none")
p2 <- p2 + theme(legend.position="none")
top <- plot_grid(p1, p2, nrow=1, align="h", axis="l", rel_widths=c(1,.9))
lgd <- plot_grid(l1, l2, ncol=1, align="v", rel_heights=c(1,4))
p <- plot_grid(top, lgd, nrow=1, align="v", rel_widths=c(1,.25))
ggsave(file.path(path, "kang_tsne.pdf"), 
    plot=p, width=14, height=5.8, unit="cm")

# ------------------------------------------------------------------------------
# nb. of cells by cluster & sample
n_cells <- data.frame(rowData(daf)) %>% group_by_("cluster_id", "sample_id") %>% tally
cids <- unique(n_cells$cluster_id)
sids <- unique(n_cells$sample_id)
n_cells <- matrix(n_cells$n,
    nrow=nlevels(sample_ids(daf)), 
    ncol=nlevels(cluster_ids(daf)),
    dimnames=list(sids, cids))
n_cells <- t(n_cells)/colSums(n_cells)*100
counts <- as.numeric(table(cluster_ids(daf)))
counts <- rowAnnotation(
    show_annotation_name=TRUE,
    annotation_name_offset=unit(1, "cm"),
    annotation_name_rot=c(0, 0, 90),
    annotation_name_gp=gpar(fontsize=8),
    "# cells"=row_anno_barplot(axis_gp=gpar(fontsize=6),
        x=counts, axis=TRUE, border=FALSE, bar_width=.8, 
        gp=gpar(fill="grey", col="white")), width=unit(1.5, "cm"))
hm <- Heatmap(n_cells,
    name="% cells\n(row normalised)",
    col=RColorBrewer::brewer.pal(9, "PRGn"),
    rect_gp=gpar(col="white"),
    column_order=levels(sample_ids(daf)),
    cluster_rows=FALSE, cluster_columns=FALSE,
    row_names_gp=gpar(fontsize=8),
    column_names_gp=gpar(fontsize=8),
    heatmap_legend_param=list(
        title_gp=gpar(fontsize=8, fontface="bold", lineheight=.8),
        labels_gp=gpar(fontsize=6)))

pdf(file.path(path, "kang_nb_cells.pdf"), width=14/2.53, height=5/2.53)
print(hm + counts)
decorate_heatmap_body("% cells\n(row normalised)", {
    grid.lines(8/16, c(0,1), gp=gpar(lwd=1.5, lty=2, col="black"))
})
dev.off()
