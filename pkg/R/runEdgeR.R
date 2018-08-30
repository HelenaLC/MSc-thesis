#' @rdname runEdgeR
#' @title ...
#'
#' @description ...
#'
#' @param x a \code{daFrame}.
#' @param k numeric or character string. Specifies the clustering to use.
#' @param method a character string.
#'
#' @author Helena Lucia Crowell
#'
#' @importFrom CATALYST cluster_codes sample_ids state_markers
#' @importFrom dplyr %>% group_by_ summarise_at
#' @importFrom edgeR DGEList estimateDisp glmQLFit glmQLFTest
#' @importFrom limma makeContrasts
#' @importFrom S4Vectors metadata
#' @importFrom scran computeSumFactors
#' @importFrom stats model.matrix
#'
#' @export

setMethod(f="runEdgeR",
    signature=signature(x="daFrame"),
    def=function(x, k=NULL, 
        method=c("rawCounts", "normCounts", "scaledCPM")) {

        md <- metadata(x)$experiment_info
        design <- model.matrix(~ 0 + md$condition)
        colnames(design) <- levels(md$condition)
        contrast <- paste(levels(md$condition), collapse="-")
        contrast <- makeContrasts(contrasts=contrast, levels=design)

        # get cluster IDs
        if (is.null(k))
            k <- names(cluster_codes(x))[1]
        k <- CATALYST:::check_validity_of_k(x, k)
        cluster_ids <- cluster_codes(x)[cluster_ids(x), k]
        
        n_clusters <- nlevels(cluster_ids)
        n_samples <- nlevels(sample_ids(x))
        df <- data.frame(
            sample_id=sample_ids(x),
            cluster_id=cluster_ids)
        
        counts <- switch(match.arg(method),
            # pseudobulks of raw counts
            rawCounts={
                cbind(df, assays(x)$counts) %>%
                    group_by_(~cluster_id, ~sample_id) %>%
                    summarise_at(state_markers(x), sum)
            },
            # pseudobulks of normalized counts
            normCounts={
                cbind(df, assays(x)$normcounts) %>%
                    group_by_(~cluster_id, ~sample_id) %>%
                    summarise_at(state_markers(x), sum)
            },
            # scaled pseudobulks of CPMs
            scaledCPM={
                # compute CPM
                cpm <- scater::calculateCPM(assays(x)$counts)
                
                # compute library-sizes by cluster & sample
                counts <- cbind(df, assays(x)$counts) %>%
                    group_by_(~cluster_id, ~sample_id) %>%
                    summarise_at(state_markers(x), sum)
                lib_sizes <- rowSums(counts[, -c(1:2)])

                # scale w/ total library size in millions
                counts <- cbind(df, cpm) %>%
                    group_by_(~cluster_id, ~sample_id) %>%
                    summarise_at(state_markers(x), sum)
                counts[, state_markers(x)] <- 
                    counts[, state_markers(x)] * lib_sizes/1e6
                counts
            }
        )
        
        # split by cluster & reformat
        counts <- split(counts, rep(seq_len(n_clusters), each=n_samples))
        counts <- lapply(counts, subset, select=-c(cluster_id, sample_id))
        counts <- lapply(counts, data.frame)
        counts <- lapply(counts, t)
        counts <- lapply(counts, `colnames<-`, levels(sample_ids(x)))

        # for ea. cluster, run DEA w/ edgeR
        n_state <- length(state_markers(x))
        res <- lapply(seq_len(n_clusters), function(k) {
            y <- edgeR::DGEList(counts[[k]])
            y <- edgeR::estimateDisp(y, design)
            fit <- edgeR::glmQLFit(y, design)
            qlf <- edgeR::glmQLFTest(fit, contrast=contrast)
            data.frame(
                cluster_id=rep(k, n_state),
                marker=state_markers(x),
                p_val=qlf$table$PValue)
        })
        res <- do.call(rbind, res)
        
        # reorder
        o <- c(sapply(seq_len(n_state), function(i) 
            seq(i, n_clusters*n_state, n_state)))
        res <- res[o, ]
        rownames(res) <- NULL
        return(res)
    }
)
