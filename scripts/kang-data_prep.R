# Kang et al., data preparation
# ==============================================================================
# This script prepares the raw data (accession # GSE96583)
# for differential analysis:
# - Seurat preprocessing
# - construction of a daFrame
# ------------------------------------------------------------------------------

# load packages
suppressPackageStartupMessages({
    library(CATALYST)
    library(DESeq2)
    library(dplyr)
    library(Matrix)
    library(Seurat)
    library(SingleCellExperiment)
    library(SummarizedExperiment)
})

# load data
message("Loading data...")
df <- read.table("kang/GSE96583_batch2_df.tsv",    sep="\t", header=TRUE)
gs <- read.table("kang/GSE96583_batch2_genes.tsv", sep="\t", stringsAsFactors=FALSE, row.names=1)
gs <- setNames(gs[, 1], rownames(gs))

x <- as.matrix(cbind(
    Matrix::readMM("kang/GSM2560248_2.1.mtx"), 
    Matrix::readMM("kang/GSM2560249_2.2.mtx")))
rownames(x) <- names(gs)
colnames(x) <- rownames(df)

# remove undetected genes
message("Removing undetected genes...")
x <- x[rowSums(x > 0) > 0, ]
gs <- gs[rownames(x)]

# preprocessing
message("Preprocessing...")
x <- CreateSeuratObject(x, min.cells=100, min.genes=200)
gs <- gs[rownames(x@raw.data)]
    
mito_gs <- grep("^MT-", rownames(x@raw.data), value=TRUE)
perc_mito <- colSums(x@raw.data[mito_gs, ]) / colSums(x@raw.data)
x <- AddMetaData(x, perc_mito, "perc_mito")
x <- FilterCells(x, subset.names=c("nGene", "perc_mito"), 
    low.thresholds=c(200,-Inf), high.thresholds=c(2500,.05))
x <- NormalizeData(x, normalization.method="LogNormalize", scale.factor=1e4, display.progress=FALSE)
x <- FindVariableGenes(x, mean.function=ExpMean, dispersion.function=LogVMR, 
    x.low.cutoff=.0125, x.high.cutoff=3, y.cutoff=.5, do.plot=FALSE, display.progress=FALSE)

# filter multiplets & unassigned cells
cells_keep <- df$multiplets == "singlet" & !is.na(df$cell)
cells_keep <- rownames(df)[cells_keep]
x <- SubsetData(x, cells.use=cells_keep, subset.raw=TRUE)
df <- df[x@cell.names, ]

# construct daFrame
message("Constructing daFrame...")
row_data <- data.frame(
    row.names=x@cell.names,
    condition=df$stim,
    patient_id=factor(df$ind),
    sample_id=paste0(df$stim, df$ind),
    cluster_id=df$cell)

col_data <- data.frame(row.names=names(gs), marker_name=names(gs), symbol=gs,
    marker_class=factor(c("type", "state")[as.numeric(rownames(x@data) %in% x@var.genes)+1], levels=c("type", "state")))

ei <- data.frame(
    condition=rep(levels(row_data$condition), 
        each=nlevels(row_data$patient_id)),
    patient_id=rep(levels(row_data$patient_id), 2))
ei$sample_id <- apply(ei, 1, paste, collapse="")

n_cells <- data.frame(row_data %>% group_by(sample_id) %>% tally)
n_cells <- setNames(n_cells[, 2], n_cells[, 1])

# subsample cells for t-SNE
set.seed(1)
tsne_inds <- sample(seq_along(x@cell.names), 2e3)
md <- list(experiment_info=ei, n_cells=n_cells,
    tsne_inds=tsne_inds, tsne=list(Y=df[tsne_inds, c("tsne1", "tsne2")]),
    cluster_codes=data.frame("20"=levels(row_data$cluster_id), check.names=FALSE))

# construct daFrame
message("Saving result...")
daf <- new("daFrame", 
    SummarizedExperiment(
        assays=list(
            counts=t(x@raw.data),
            normcounts=as.matrix(t(exp(1)^x@data-1)), 
            exprs=as.matrix(t(x@data))),
        rowData=row_data, colData=col_data, metadata=md))
saveRDS(daf, "results/kang/kang_daf.rds")
