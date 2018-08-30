#' @rdname calcMeans
#' @title Calculate cluster means
#' 
#' @description Calculate cluster means 
#' (mean expression for each cluster-sample-marker combination).
#' 
#' @details Calculate mean marker expression for each cluster and sample 
#' (i.e. means for each cluster-sample-marker combination).
#' 
#' The data object is assumed to contain a factor \code{marker_class} in the column
#' meta-data (see \code{\link{prepareData}}), which indicates the protein marker class for
#' each column of data (\code{"type"}, \code{"state"}, or \code{"none"}).
#' 
#' The cluster means are required for testing for differential states within cell
#' populations, and for plotting purposes.
#' 
#' Variables \code{id_type_markers} and \code{id_state_markers} are saved in the
#' \code{metadata} slot of the output object. These can be used to identify the 'cell
#' type' and 'cell state' markers in the list of \code{assays} in the output
#' \code{\link{SummarizedExperiment}} object, which is useful in later steps of the
#' 'diffcyt' pipeline.
#' 
#' Results are returned as a new \code{\link{SummarizedExperiment}} object, where rows =
#' clusters, columns = samples, sheets (\code{assays} slot) = markers. Note that there is
#' a separate table of values (\code{assay}) for each marker. The \code{metadata} slot
#' also contains variables \code{id_type_markers} and \code{id_state_markers}, which can
#' be used to identify the sets of cell type and cell state markers in the list of
#' \code{assays}.
#' 
#' 
#' @param d_se Data object from previous steps, in \code{\link{SummarizedExperiment}}
#'   format, containing cluster labels as a column in the row meta-data (from
#'   \code{\link{generateClusters}}). Column meta-data is assumed to contain a factor
#'   \code{marker_class}.
#' 
#' 
#' @return \code{d_means}: \code{\link{SummarizedExperiment}} object, where rows =
#'   clusters, columns = samples, sheets (\code{assays} slot) = markers. The
#'   \code{metadata} slot contains variables \code{id_type_markers} and
#'   \code{id_state_markers}, which can be accessed with
#'   \code{metadata(d_means)$id_type_markers} and
#'   \code{metadata(d_means)$id_state_markers}.
#' 
#' 
#' @importFrom SummarizedExperiment SummarizedExperiment assay rowData colData
#' @importFrom dplyr group_by_ tally summarize
#' @importFrom tidyr complete
#' @importFrom reshape2 acast
#' @importFrom magrittr '%>%'
#' @importFrom methods is
#' 
#' @export

setMethod(f="calcMeans",
    signature=signature(d_se="SummarizedExperiment"),
    definition=function(d_se) {
        
        if (!("cluster_id" %in% (colnames(rowData(d_se)))))
            stop("Data object does not contain cluster labels.\n", 
                " Run 'generateClusters' to generate cluster labels.")
        
        is_marker <- colData(d_se)$marker_class != "none"
        
        # identify 'cell type' and 'cell state' markers in final list of assays
        id_type_markers <- (colData(d_se)$marker_class == "type")[is_marker]
        id_state_markers <- (colData(d_se)$marker_class == "state")[is_marker]
        
        # calculate cluster means for each marker
        
        assaydata_mx <- assay(d_se)
        
        means <- vector("list", sum(is_marker))
        marker_names_sub <- as.character(colData(d_se)$marker_name[is_marker])
        names(means) <- marker_names_sub
        
        clus <- rowData(d_se)$cluster_id
        smp <- rowData(d_se)$sample_id
        
        for (i in seq_along(means)) {
            assaydata_i <- assaydata_mx[, marker_names_sub[i], drop = FALSE]
            assaydata_i <- as.data.frame(assaydata_i)
            assaydata_i <- cbind(assaydata_i, sample_id = smp, cluster_id = clus)
            colnames(assaydata_i)[1] <- "value"
            
            assaydata_i %>% 
                group_by_(~cluster_id, ~sample_id) %>% 
                summarize(mean = mean(value)) -> 
                mean
            
            mean <- reshape2::acast(mean, cluster_id ~ sample_id, value.var = "mean", fill = NA)
            
            # fill in any missing clusters
            if (nrow(mean) < nlevels(rowData(d_se)$cluster_id)) {
                ix_missing <- which(!(levels(rowData(d_se)$cluster_id) %in% rownames(mean)))
                mean_tmp <- matrix(NA, nrow = length(ix_missing), ncol = ncol(mean))
                rownames(mean_tmp) <- ix_missing
                mean <- rbind(mean, mean_tmp)
                # re-order rows
                mean <- mean[order(as.numeric(rownames(mean))), , drop = FALSE]
            }
            
            means[[i]] <- mean
        }
        
        # check cluster IDs and sample IDs are identical
        for (i in seq_along(means)) {
            if (!all(rownames(means[[i]]) == rownames(means[[1]]))) {
                stop("Cluster IDs do not match")
            }
            if (!all(colnames(means[[i]]) == colnames(means[[1]]))) {
                stop("Sample IDs do not match")
            }
        }
        
        # create new SummarizedExperiment (rows = clusters, columns = samples)
        
        row_data <- data.frame(
            cluster_id = factor(rownames(means[[1]]), levels = levels(rowData(d_se)$cluster_id)), 
            stringsAsFactors = FALSE
        )
        
        col_data <- metadata(d_se)$experiment_info
        
        # rearrange sample order to match 'experiment_info'
        means <- lapply(means, function(m) {
            m[, match(col_data$sample_id, colnames(m)), drop = FALSE]
        })
        stopifnot(all(sapply(means, function(m) {
            col_data$sample_id == colnames(m)
        })))
        
        metadata <- list(id_type_markers = id_type_markers, 
            id_state_markers = id_state_markers)
        
        d_means <- SummarizedExperiment(
            assays = list(means = means), 
            rowData = row_data, 
            colData = col_data, 
            metadata = metadata
        )
        
        d_means
    }
)
        