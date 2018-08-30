# simDD, simulation example
# ==============================================================================
# Generate figures for an examplary scDD simulation:
# - expression distributions
# - expression heatmap
# ------------------------------------------------------------------------------

path <- "results/simDD-sim_ex"

# load packages
suppressPackageStartupMessages({
    library(CATALYST)
    library(ComplexHeatmap)
    library(dplyr)
    library(ggplot2)
    library(methods)
    library(pkg)
    library(RColorBrewer)
    library(SummarizedExperiment)
})

# load data
data("sceEx", package="pkg")
data("idxEx", package="pkg")

# simulate data
sim <- simDD2(sceEx, idxEx, 
    n_cells=6e3, n_genes=1e3, n_clusters=1, n_samples=3,
    perc_type=0.2, perc_de=1, cat_probs=rep(0.25, 4), seed=123)

# expression profiles
# ------------------------------------------------------------------------------
es <- log2(assays(sim)$counts+1)
colnames(es) <- paste0(colnames(es), "(", metadata(sim)$gene_info$cat, ")")

# sample 4 type genes + 8 state genes (2 per category)
set.seed(1)
genes_keep <- c(
    sample(which(colData(sim)$marker_class == "type" & metadata(sim)$gene_info$cat == "EE"), 2),
    sample(which(colData(sim)$marker_class == "type" & metadata(sim)$gene_info$cat == "EP"), 2),
    replicate(2, c(
        sample(which(metadata(sim)$gene_info$cat == "DE"), 1),
        sample(which(metadata(sim)$gene_info$cat == "DP"), 1), 
        sample(which(metadata(sim)$gene_info$cat == "DM"), 1),
        sample(which(metadata(sim)$gene_info$cat == "DB"), 1))))

# plot results
df <- data.frame(es[, genes_keep], rowData(sim), check.names=FALSE)
df <- subset(df, select=- cluster_id)
df <- reshape2::melt(df, value.name="expression")
p <- ggplot(df, aes(x=expression, group=sample_id, col=condition)) +
    facet_wrap(~ variable, scales="free", ncol=4) +
    scale_x_continuous(expand=c(.05,0)) + 
    scale_y_continuous(expand=c(.075,0)) +
    geom_density(size=0.5) + pkg:::prettify + 
    theme(aspect.ratio=2/3,
        strip.text=element_text(hjust=0.5),
        legend.position="none")
ggsave(plot=p, width=14, height=9, unit="cm", 
    file.path(path, "simDD_ex-expr_profiles.pdf"))

# heatmap of median expressions
# ------------------------------------------------------------------------------

# simulate data
sim <- simDD(x=sceEx, y=idxEx, 
    n_cells=3e3, n_genes=300, n_clusters=6, n_samples=3,
    perc_type=0.5, perc_de=0.5, cat_probs=rep(0.25, 4), seed=1)

# samples genes for plotting
set.seed(88)
cluster_id <- 1
gi <- metadata(sim)$gene_info
gi <- gi[gi$cluster_id == cluster_id, ]
gi$gene <- as.character(gi$gene)
t_gs <- gsub("[0-9]", "", gi$gene) == "type"
s_gs <- gsub("[0-9]", "", gi$gene) == "state"
n <- 4
genes_keep <- c(
    sample(gi$gene[t_gs & gi$cat == "EE"], n),
    sample(gi$gene[t_gs & gi$cat == "EP"], n),
    sapply(levels(gi$cat), function(cat)
        sample(gi$gene[s_gs & gi$cat == cat], n)))
gi <- gi[gi$gene %in% genes_keep, ]

# plot median exrps. by sample for exemplary cluster
cells_keep <- cluster_ids(sim) == cluster_id
es <- log2(assays(sim)$counts+1)
es <- es[cells_keep, genes_keep]
es <- CATALYST:::scale_exprs(es)
df <- data.frame(es, sample_id=sample_ids(sim)[cells_keep])
df <- df %>% group_by_("sample_id") %>% summarise_all(median)
df <- t(data.frame(df, row.names=1))
gi <- gi[match(rownames(df), gi$gene), ]

hm <- Heatmap(
    mat=df, name="median scaled\nexpression",
    col=rev(brewer.pal(11, "RdYlBu")),
    rect_gp=gpar(col="white"),
    column_order=levels(sample_ids(sim)),
    cluster_rows=FALSE, cluster_columns=FALSE,
    column_title="sample_id", column_title_side="bottom",
    split=gi$cat, gap=unit(.5, "cm"),
    column_title_gp=gpar(fontsize=8, fontface="bold"), 
    row_names_gp=gpar(fontsize=6, rot=0),
    row_title_gp=gpar(fontsize=8),
    column_names_gp=gpar(fontsize=8),
    heatmap_legend_param=list(at=seq(0, 1, 0.2),
        title_gp=gpar(fontsize=8, fontface="bold", lineheight=0.8),
        labels_gp=gpar(fontsize=6)))
row_anno <- rowAnnotation(
    df=data.frame(is_de=gi$is_de, row.names=gi$gene),
    col=list(is_de=c("TRUE"="limegreen", "FALSE"="grey")),
    annotation_legend_param=list(
        title_gp=gpar(fontsize=8, fontface="bold"),
        labels_gp=gpar(fontsize=6)),
    gp=gpar(col="white"))
pdf(sprintf(file.path(path, "simDD_ex-med_exprs.pdf"), cluster_id),
    width=14/2.53, height=10/2.53)
hm + row_anno   
dev.off()
