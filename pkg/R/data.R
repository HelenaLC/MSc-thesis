# ==============================================================================
#' @rdname data
#' @name data
#' @aliases idxEx sceEx
#' @title Example data set
#' 
#' @description 
#' \code{sceEx}\describe{
#' Toy example \code{SingleCellExperiment} from the \code{scDD} package with
#' 500 genes and 142 samples (78 from condition 1 & 64 from condition 2).
#' Condition labels (1 or 2) are stored in the colData slot.}
#' 
#' \code{idxEx}\describe{
#' Vector of indices of a reasonable set of genes to use for the simulation.}
#' 
#' @return see descriptions above.
#' 
#' @examples
#' data(sceEx) # SingleCellExperiment
#' data(idxEx) # gene indices
#' 
#' @references
#' Korthauer KD, Chu LF, Newton MA, Li Y, Thomson J, Stewart R, Kendziorski C. 
#' A statistical approach for identifying differential distributions in single-
#' cell RNA-seq experiments. \emph{Genome Biology}. 2016 Oct 25;17(1):222. 
#' https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1077-y
#' 
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
# ------------------------------------------------------------------------------
NULL