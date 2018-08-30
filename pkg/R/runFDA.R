#' @rdname runFDA
#' @title ...
#'
#' @description ...
#'
#' @param x a \code{daFrame}.
#' @param k numeric or character string. Specifies the clustering to use.
#' @param reso resolution.
#' @param n_perm numeric. Nb. of permutations to use in creating
#'   the null distribution. Passed to \code{\link[fda]{tperm.fd}}.
#'
#' @author Helena Lucia Crowell
#'
#' @importFrom CATALYST cluster_codes sample_ids state_markers
#' @importFrom data.table data.table
#' @importFrom fda Data2fd tperm.fd
#' @importFrom gtools permutations
#' @importFrom methods is
#' @importFrom stats ecdf
#'
#' @export

setMethod(f="runFDA",
    signature=signature(x="daFrame"),
    def=function(x, k=NULL, reso=50, n_perm=200) {
        
        # get cluster IDs
        if (is.null(k))
            k <- names(cluster_codes(x))[1]
        k <- CATALYST:::check_validity_of_k(x, k)
        cluster_ids <- cluster_codes(x)[cluster_ids(x), k]
        cluster_nms <- levels(cluster_ids)
        
        # split exprs. by cluster & sample
        es <- log2(assays(x)$counts+1)[, state_markers(x)]
        dt <- data.table(es, cluster_id=cluster_ids)
        es_by_cluster <- split(dt, by="cluster_id", keep.by=FALSE)
        es_by_cluster_sample <- lapply(cluster_nms, function(k) {
            es <- es_by_cluster[[k]]
            sample_ids <- sample_ids(x)[cluster_ids == k]
            dt <- data.table(es, sample_id=sample_ids)
            split(dt, by="sample_id", keep.by=FALSE)
        })
        names(es_by_cluster_sample) <- cluster_nms

        # get expression range by cluster & gene
        rs <- lapply(es_by_cluster, apply, 2, range, finite = TRUE)
        xs <- lapply(rs, apply, 2, function(r) seq(r[1], r[2], length=reso))

        # get ECDFs by cluster, sample & gene
        sample_ids <- levels(sample_ids(x))
        
        ecdfs <- lapply(es_by_cluster_sample, sapply, sapply, stats::ecdf)
        ecdfs.x <- lapply(cluster_nms, function(k) lapply(sample_ids, function(s)
            sapply(state_markers(x), function(g) ecdfs[[k]][, s][[1]](xs[[k]][, g]))))
        names(ecdfs.x) <- cluster_nms

        md <- metadata(x)$experiment_info
        groups <- levels(md$condition)
        g1 <- grepl(groups[1], sample_ids)

        combis <- expand.grid(stringsAsFactors=FALSE,
            cluster_id=cluster_nms, marker=state_markers(x))
        
        n_samples <- length(sample_ids)
        n_perms <- ncol(combn(n_samples, n_samples/2))
        if (n_perm > n_perms) {
            message("Parameter 'n_perm' exceeds the number of distinct",
                "  permutations. Reducing 'n_perm' to ", paste0(n_perms, "."))
            n_perm <- n_perms
        }
        
        pvals <- mapply(function(k, g) {
            ecdfs_g1 <- sapply(ecdfs.x[[k]][ g1], function(x) x[, g])
            ecdfs_g2 <- sapply(ecdfs.x[[k]][!g1], function(x) x[, g])
            
            fd1 <- fda::Data2fd(argvals=xs[[k]][, g], y=ecdfs_g1)
            fd2 <- fda::Data2fd(argvals=xs[[k]][, g], y=ecdfs_g2)
            tp <- tryCatch({
                fda::tperm.fd(fd1, fd2, nperm=n_perm,
                    argvals=xs[[k]][, g], plotres=FALSE)
            }, error=function(e) return(list(error=e)))
            err <- is(tp$error, "simpleError")
            if (err) return(NA)
            return(tp$pval)
        }, k=combis$cluster_id, g=combis$marker)
        
        pvals <- pvals + 1 / (n_perm + 1)
        pvals[pvals > 1] <- 1
        pvals <- data.frame(combis, p_val=pvals)
        return(pvals)
    }
)

