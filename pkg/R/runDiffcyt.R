#' @rdname runDiffcyt
#' @title ...
#'
#' @description ...
#'
#' @param x a \code{daFrame}.
#' @param k numeric or character string.
#'   Specifies the clustering from which to retrieve cluster IDs.
#'   Must be one of \code{colnames(cluster_codes(x))}.
#'   Defaults to the 1st clustering available.
#' @param method a character string.
#' @param data character string specifying the data to use.
#'   Should be one of \code{assayNames(x)}.
#' @param log logical. Should log2 be taken before computing means/medians.
#' @param fun character string specifying the summary statistic to use.
#'
#' @author Helena Lucia Crowell
#'
#' @import SingleCellExperiment
#' @importFrom CATALYST cluster_codes cluster_ids sample_ids
#' @importFrom diffcyt calcCounts calcMedians testDS_limma testDS_LMM
#'   createContrast createDesignMatrix createFormula
#' @importFrom dplyr %>% group_by summarize
#' @importFrom S4Vectors metadata
#' @importFrom scater calculateCPM
#' @importFrom scran computeSumFactors
#' @importFrom SummarizedExperiment assayNames
#'
#' @export

setMethod(f = "runDiffcyt",
    signature = signature(x = "daFrame"),
    definition = function(x, data, log = FALSE, k = NULL, 
        method = c("limma", "LMM"), fun = c("mean", "median")) {
        
        method <- match.arg(method)
        fun <- match.arg(fun)
        
        cs <- diffcyt::calcCounts(x)
        stopifnot(data %in% assayNames(x))
        assay(x) <- assays(x)[[data]]
        if (log) 
            assay(x) <- log2(assay(x) + 1)
        
        ms <- switch (fun,
            mean = pkg::calcMeans(x),
            median = diffcyt::calcMedians(x))
        
        contrast <- createContrast(c(0, 1))
        md <- metadata(x)$experiment_info
        col <- grep("condition", colnames(md))
        
        res <- switch(method,
            limma = {
                design <- createDesignMatrix(md, cols_design = col)
                testDS_limma(cs, ms, design, contrast)
            },
            LMM = {
                formula <- createFormula(md, cols_fixed = col)
                testDS_LMM(cs, ms, formula, contrast)
            })
        rowData(res)
    }
)
