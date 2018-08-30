#' @rdname calcMeans
#' @param ... optional arguments.
#' @export
setGeneric("calcMeans",
    function(d_se) standardGeneric("calcMeans"))

#' @rdname runDiffcyt
#' @param ... optional arguments.
#' @export
setGeneric("runDiffcyt",
    function(x, ...) standardGeneric("runDiffcyt"))

#' @rdname runEdgeR
#' @param ... optional arguments.
#' @export
setGeneric("runEdgeR",
    function(x, ...) standardGeneric("runEdgeR"))

#' @rdname runFDA
#' @param ... optional arguments.
#' @export
setGeneric("runFDA",
    function(x, ...) standardGeneric("runFDA"))

#' @rdname simDD
#' @param ... optional arguments.
#' @export
setGeneric("simDD",
    function(x, y, ...) standardGeneric("simDD"))

#' @rdname simDD2
#' @param ... optional arguments.
#' @export
setGeneric("simDD2",
    function(x, y, ...) standardGeneric("simDD2"))

#' @rdname plotPerf
#' @param ... optional arguments.
#' @export
setGeneric("plotPerf",
    function(x, ...) standardGeneric("plotPerf"))