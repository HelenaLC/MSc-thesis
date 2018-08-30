#' @rdname simDD
#' @title \code{scDD} simulation
#'
#' @description Complex \code{scDD} simulation with 
#' varying numbers and sizes of samples and clusters.
#'
#' @param x a \code{SingleCellExperiment}.
#' @param y numeric vector of gene indices to simulate from.
#' @param n_cells numeric of length 1, 2, or \code{n_clusters}.
#'   Specifies the nb. of cells per cluster. 
#'   If a single value is provided, clusters will be equally-sized.
#'   If a range is provided and \code{rand_clusters = TRUE}, 
#'   cluster sizes will be sampled from therein.
#'   If \code{n_clusters} values are provided, these will be used.
#' @param n_genes numeric. Nb. of genes.
#' @param n_samples numeric. Nb. of samples per condition.
#' @param n_clusters numeric. Nb. of clusters.
#' @param rand_samples,rand_clusters logical. 
#'   Should sample/cluster sized be randomized?
#' @param perc_type numeric. Percentage of genes to set as type genes.
#' @param perc_de numeric of length 1 or 2. 
#'   Specifies the percentage of state genes to be DE per cluster.
#'   If a range is provided, the percentage will be sampled from therein.
#' @param cat_probs numeric of length 4 specifying the probabilities 
#'   of DE genes to be of category DE, DP, DM, and DB, respectively.
#' @param sd_range numeric vector of length to specifying
#'   the range of SDs of FCs to sample from.
#' @param mode_fc numeric vector of possible FCs b/w modes for DP, DM, and DB.
#' @param zeros method to generate zeros. 
#' @param seed numeric. Random seed for sampling & \code{scDD} simulation.
#'
#' @return a \code{daFrame}.
#' 
#' @author Helena Lucia Crowell
#'
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @importFrom BiocGenerics counts
#' @importFrom CATALYST daFrame
#' @importFrom methods new
#' @importFrom scater calculateCPM calculateQCMetrics isOutlier normalise
#' @importFrom scDD simulateSet
#' @importFrom stats runif setNames
#' @importFrom zeallot %<-%
#'
#' @examples
#' sce <- get(data("scDatEx", package = "scDD"))
#' gs <- scDD:::findIndex(sce, condition = "condition")
#' simDD(x = sce, y = gs)
#'
#' @export

setMethod(f="simDD",
    signature=signature(x="SingleCellExperiment", y="numeric"),
    def=function(x, y, n_cells=100, n_genes=5, 
        n_samples=4, n_clusters=1, rand_samples=FALSE, rand_clusters=FALSE, 
        perc_type=0, perc_de=c(0, .25), cat_probs=rep(.25, 4), 
        sd_range=c(1, 3), mode_fc=c(2, 3, 4), seed=1) {
        
        # check validity of input arguments
        args <- as.list(environment())
        checkArgs_simDD(args)
        
        # get nb. of cells per cluster
        if (length(n_cells) == 1) {
            cells_per_cluster <- rep(n_cells, n_clusters)
        } else if (length(n_cells) == 2 & rand_clusters) {
            set.seed(seed)
            cells_per_cluster <- sample(n_cells[1]:n_cells[2], n_clusters)
        } else if (length(n_cells) == n_clusters) {
            cells_per_cluster <- n_cells
        }
        
        is_odd <- sapply(cells_per_cluster, function(n) n %% 2 != 0)
        cells_per_cluster[is_odd] <- cells_per_cluster[is_odd]+1
        
        # get nb. of type & state genes
        n_type <- round(n_genes * perc_type)
        n_state <- n_genes - n_type
        
        message("Simulating...\n", 
            " * ", n_genes, " genes", 
            " (", n_type, " type, ", n_state, " state)\n",
            " * ", n_clusters, " clusters",
            " (", paste(cells_per_cluster, collapse=", "), " cells per cluster)")
        
        # initialise count & gene info. matrix
        type_gs <- sprintf("type%s", seq_len(n_type))
        state_gs <- sprintf("state%s", seq_len(n_state))
        gs <- c(type_gs, state_gs)
        counts <- matrix(NA, n_genes, sum(cells_per_cluster), dimnames=list(gs, NULL))
        
        de_cats <- c("DE", "DP", "DM", "DB")
        de_cats <- factor(de_cats, levels=de_cats)
        
        cats <- c("DE", "DP", "DM", "DB", "EE", "EP")
        cats <- factor(character(n_genes*n_clusters), levels=cats)
        fcs <- numeric(n_genes*n_clusters)
        
        # simulate n_clusters
        inds <- c(0, cumsum(cells_per_cluster))
        inds <- lapply(seq_along(inds)[-1], function(k) (inds[k-1]+1):inds[k])
        for (k in seq_len(n_clusters)) {
            set.seed(seed+(k-1))
            
            # sample nb. of DE genes
            p_de <- ifelse(length(perc_de) == 1, 
                perc_de, runif(1, perc_de[1], perc_de[2]))
            n_de <- round(n_state * p_de)
            n_non_de <- n_genes - n_de
            
            # sample gene categories
            n_ee <- sample(n_non_de, 1)
            n_ep <- n_non_de - n_ee
            n_de <- c(table(sample(de_cats, n_de, TRUE, cat_probs)))
            
            # simulate cluster k
            sim <- scDD::simulateSet(
                sce=x, gs=y, numSamples=cells_per_cluster[k]/2,
                n_de[1], n_de[2], n_de[3], n_de[4], n_ee, n_ep, 
                sd.range=sd_range, modeFC=mode_fc, seed=seed+(k-1))
            
            # sample type & state genes
            set.seed(seed+(k-1))
            idx_type <- which(!rowData(sim)$Category %in% de_cats)
            idx_type <- sample(idx_type, n_type)
            idx_state <- sample(setdiff(seq_len(n_genes), idx_type))
            
            counts[type_gs, inds[[k]]] <- assay(sim)[idx_type, ]
            counts[state_gs, inds[[k]]] <- assay(sim)[idx_state, ]
            
            # store gene categories & FCs
            idx <- seq_len(n_genes)+(k-1)*n_genes
            cats[idx] <- rowData(sim)$Category[c(idx_type, idx_state)]
            fcs[idx] <- rowData(sim)$FC[c(idx_type, idx_state)]
        }
        gi <- data.frame(cluster_id=rep(seq_len(n_clusters), each=n_genes),
            gene=rep(gs, n_clusters), cat=cats, fc=fcs, is_de=cats %in% de_cats)
        
        genes_keep <- apply(counts, 1, function(g) !any(is.na(g)))
        gi <- gi[rep(genes_keep, n_clusters), ]
        counts <- counts[genes_keep, ]
        
        # normalise
        message("Normalising...")
        cpm <- scater::calculateCPM(counts, use_size_factors=FALSE)
        
        # get nb. of cells per sample
        if (rand_samples) {
            set.seed(seed)
            cells_per_sample <- 
                replicate(2, lapply(seq_len(n_clusters), function(k) {
                    n <- cells_per_cluster[k]/2
                    mins <- floor(n/(n_samples+1))
                    maxs <- ceiling(n/n_samples)
                    ns <- replicate(n_samples-1, sample(mins:maxs, 1))
                    c(ns, n - sum(ns))
                }))
        } else {
            cells_per_sample <- 
                replicate(2, lapply(seq_len(n_clusters), function(k) {
                    n <- cells_per_cluster[k]/2
                    ns <- rep(floor(n/n_samples), n_samples-1)
                    c(ns, n - sum(ns))
                }))
        }
        if (n_clusters == 1) {
            cells_per_sample <- matrix(unlist(cells_per_sample))
        } else {
            cells_per_sample <- apply(cells_per_sample, 1, unlist)
        }
        
        # rowData
        condition <- unlist(lapply(cells_per_cluster, 
            function(n) rep(c(1, 2), each=n/2)))
        condition <- c("ref", "trt")[condition]
        sample_id <- unlist(apply(cells_per_sample, 2, 
            function(n) rep.int(rep(seq_len(n_samples), 2), n)))
        sample_id <- paste0(condition, sample_id)
        cluster_id <- rep.int(seq_len(n_clusters), cells_per_cluster)
        cluster_id <- factor(cluster_id, levels=seq_len(n_clusters))
        row_data <- data.frame(sample_id, condition, cluster_id)
        
        # colData
        col_data <- data.frame(marker_name=gs, 
            marker_class=factor(rep.int(c("type", "state"), c(n_type, n_state))))
        col_data <- col_data[genes_keep, ]
        
        # metadata
        condition <- rep(c("ref", "trt"), each=n_samples)
        sample_id <- paste0(condition, rep(seq_len(n_samples), 2))
        md <- list(
            experiment_info=data.frame(condition, sample_id),
            n_cells=table(row_data$sample_id),
            cluster_codes=data.frame("20"=factor(seq_len(n_clusters)), check.names=FALSE),
            gene_info=gi)
        
        # construct daFrame
        new("daFrame", SummarizedExperiment(
            assays=list(counts=t(counts), normcounts=t(cpm)),
            rowData=row_data, colData=col_data, metadata=md))
    }
)

# check validity of input arguments
checkArgs_simDD <- function(args) {
    for (i in seq_along(args)) assign(names(args)[i], args[[i]])
    stopifnot(is.numeric(y), min(y) > 0, max(y) <= nrow(x),
        is.numeric(n_cells), min(n_cells) > 0,
        length(n_cells) %in% c(1, 2, n_clusters),
        is.numeric(n_genes), length(n_genes) %in% c(1, 6),
        is.numeric(n_samples), length(n_samples) == 1,
        is.numeric(n_clusters), length(n_clusters) == 1, n_clusters > 0,
        is.numeric(sd_range), length(sd_range) == 2,
        is.numeric(mode_fc), !any(mode_fc <= 0),
        is.numeric(perc_type), length(perc_type) == 1,
        perc_type >= 0, perc_type <= 100,
        is.numeric(seed), length(seed) == 1)
}