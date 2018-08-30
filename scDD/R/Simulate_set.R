#' simulateSet
#'
#' Simulation of a complete dataset, where the number of each type of 
#' differential distributions and equivalent distributions is specified.
#'
#' @inheritParams singleCellSimu
#' @inheritParams scDD
#' @inheritParams findFC
#' @inheritParams findIndex
#' 
#' @param nDE Number of DE genes to simulate
#' @param nDP Number of DP genes to simulate
#' @param nDM Number of DM genes to simulate
#' @param nDB Number of DB genes to simulate
#' @param nEE Number of EE genes to simulate
#' @param nEP Number of EP genes to simulate
#' @param seed numeric.
#'  
#' @param param a \code{MulticoreParam} or \code{SnowParam} object of 
#' the \code{BiocParallel} package that defines a parallel backend.  
#' The default option is \code{BiocParallel::bpparam()} which will 
#' automatically creates a cluster appropriate for the operating system.  
#' Alternatively, the user can specify the number of cores they wish to use by
#' first creating the corresponding \code{MulticoreParam} (for Linux-like OS) 
#' or \code{SnowParam} (for Windows) object, and then passing it into the 
#' \code{scDD} function. This could be done to specify a parallel backend on 
#' a Linux-like OS with, say 12 cores by setting 
#' \code{param=BiocParallel::MulticoreParam(workers=12)}.
#' 
#' @param zeros method to generate zeros.
#' @param perc_zeros numeric between 0 and 1. 
#'   Specifies the fixed proportion of zeros in every gene.
#'   Ignored if \code{zeros != "constant"}.
#' 
#' @import SingleCellExperiment
#' @importFrom BiocParallel register
#' 
#' @export
#'
#' @references Korthauer KD, Chu LF, Newton MA, Li Y, Thomson J, Stewart R, 
#' Kendziorski C. A statistical approach for identifying differential 
#' distributions
#' in single-cell RNA-seq experiments. Genome Biology. 2016 Oct 25;17(1):222. 
#' \url{https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-
#' 1077-y}
#'  
#' @return  An object of class \code{SingleCellExperiment} that contains 
#' simulated single-cell expression and metadata. The \code{assays} 
#'   slot contains a named list of matrices, where the simulated counts are 
#'   housed in the one named \code{normcounts}.  This matrix should have one
#'    row for each gene (\code{nDE + nDP + nDM + nDB + nEE
#'    + nEP} rows) and one sample for each column (\code{numSamples} columns).  
#'   The \code{colData} slot contains a data.frame with one row per 
#'   sample and a column that represents biological condition, which is 
#'   in the form of numeric values (either 1 or 2) that indicates which 
#'   condition each sample belongs to (in the same order as the columns of 
#'   \code{normcounts}). The \code{rowData} slot contains information about the
#'   category of the gene (EE, EP, DE, DM, DP, or DB), as well as the simulated
#'   foldchange value.
#' 
#' @examples 
#' 
#' # Load toy example ExpressionSet to simulate from
#' data(scDatEx)
#' 
#' # check that this object is a member of the ExpressionSet class
#' # and that it contains 142 samples and 500 genes
#' 
#' class(scDatEx)
#' show(scDatEx)
#' 
#' # set arguments to pass to simulateSet function
#' # we will simuate 30 genes total; 5 genes of each type;
#' # and 100 samples in each of two conditions
#' 
#' nDE <- 5
#' nDP <- 5
#' nDM <- 5
#' nDB <- 5
#' nEE <- 5
#' nEP <- 5
#' numSamples <- 100
#' seed <- 816
#' 
#' # create simulated set with specified numbers of 
#' # DE, DP, DM, DM, EE, and EP genes,
#' # specified number of samples, DE genes are 2 standard deviations apart, 
#' # and multimodal genes have modal distance of 4 standard deviations
#' 
#' SD <- simulateSet(scDatEx, numSamples=numSamples, 
#'     nDE=nDE, nDP=nDP, nDM=nDM, nDB=nDB, nEE=nEE, nEP=nEP, 
#'     sd.range=c(2,2), modeFC=4, seed=seed)

simulateSet <- function(sce, gs=NULL, numSamples=100, 
    nDE=250, nDP=250, nDM=250, nDB=250, nEE=5000, nEP=4000, 
    sd.range=c(1,3), modeFC=c(2,3,4), var_inf=NULL, 
    condition="condition", seed=284, param=bpparam(),
    zeros=c("empirical", "simulated", "constant"),
    perc_zeros=0.5) {
    
    # reference category/condition - the first listed one
    ref <- unique(colData(sce)[[condition]])[1]  
    
    BiocParallel::register(BPPARAM=param)
    
    if (is.null(gs))
        gs <- findIndex(sce, condition)

    FC <- scDD:::findFC(sce, gs, sd.range=sd.range, N=6, 
        overExpressionProb=0.5, plot.FC=FALSE, condition)
    
    generateZero <- match.arg(zeros)
    constantZero <- perc_zeros
    
    # pull off matrix of expression values for condition 1
    Dataset1 <- normcounts(sce[,colData(sce)[[condition]]==ref])
    
    set.seed(seed)
    
    # initialize pe_mat and fcs
    pe_mat <- NULL
    fcs <- NULL
    
    ### DE
    if (nDE > 0){
        SD1 <- scDD:::singleCellSimu(Dataset1, Method="DE", gs, FC, modeFC, 
            Validation=FALSE, numGenes=nDE, numDE=nDE,
            numSamples=numSamples, generateZero=generateZero,
            constantZero=constantZero, varInflation=var_inf)
        Simulated_Data <- SD1[[1]]
        rnms <- rep("EE", nrow(Simulated_Data))
        rnms[SD1[[2]]] <- "DE"
        rownames(Simulated_Data) <- rnms
        pe_mat <- Simulated_Data
        fcs <- SD1[[3]]
    }
    
    ### DP
    if (nDP > 0){
        SD2 <- singleCellSimu(Dataset1, Method="DP", gs, FC, modeFC, 
            DP=c(0.33,0.66), Validation=FALSE,
            numGenes=nDP, numDE=nDP, numSamples=numSamples,
            generateZero=generateZero, constantZero=constantZero, 
            varInflation=var_inf)
        Simulated_Data_DP <- SD2[[1]]
        rnms <- rep("EE", nrow(Simulated_Data_DP))
        rnms[SD2[[2]]] <- "DP"
        rownames(Simulated_Data_DP) <- rnms
        pe_mat <- rbind(pe_mat, Simulated_Data_DP)
        fcs <- c(fcs, SD2[[3]])
    }
    
    ### DM
    if (nDM > 0){
        SD3 <- singleCellSimu(Dataset1, Method="DM", gs, FC, modeFC, 
            Validation=FALSE,
            numGenes=nDM, numDE=nDM, numSamples=numSamples, 
            generateZero=generateZero, constantZero=constantZero, 
            varInflation=var_inf)
        Simulated_Data_DM <- SD3[[1]]
        rnms <- rep("EE", nrow(Simulated_Data_DM))
        rnms[SD3[[2]]] <- "DM"
        rownames(Simulated_Data_DM) <- rnms
        pe_mat <- rbind(pe_mat, Simulated_Data_DM)
        fcs <- c(fcs, SD3[[3]])
    }
    
    ### DB
    if (nDB > 0){
        SD4 <- singleCellSimu(Dataset1, Method="DB", gs, FC, modeFC, 
            DP=c(0.5,0.5), Validation=FALSE,
            numGenes=nDB, numDE=nDB, numSamples=numSamples, 
            generateZero=generateZero, constantZero=constantZero, 
            varInflation=var_inf)
        Simulated_Data_DB <- SD4[[1]]
        rnms <- rep("EE", nrow(Simulated_Data_DB))
        rnms[SD4[[2]]] <- "DB"
        rownames(Simulated_Data_DB) <- rnms
        pe_mat <- rbind(pe_mat, Simulated_Data_DB)
        fcs <- c(fcs, SD4[[3]])
    }
    
    ### EP
    if (nEP > 0){
        SD5 <- singleCellSimu(Dataset1, Method="DP", gs, FC, modeFC, 
            DP=c(0.50,0.50), Validation=FALSE,
            numGenes=nEP, numDE=0, numSamples=numSamples, 
            generateZero=generateZero, constantZero=constantZero, 
            varInflation=var_inf)
        Simulated_Data_EP <- SD5[[1]]
        rnms <- rep("EP", nrow(Simulated_Data_EP))
        rnms[SD5[[2]]] <- "DP"
        rownames(Simulated_Data_EP) <- rnms
        pe_mat <- rbind(pe_mat, Simulated_Data_EP)
        fcs <- c(fcs, SD5[[3]])
    }
    
    ### EE
    if (nEE > 0){
        SD6<- singleCellSimu(Dataset1, Method="DE", gs, FC, modeFC, 
            Validation=FALSE, 
            numGenes=nEE, numDE=0, numSamples=numSamples, 
            generateZero=generateZero, constantZero=constantZero, 
            varInflation=var_inf)
        Simulated_Data_EE <- SD6[[1]]
        rnms <- rep("EE", nrow(Simulated_Data_EE))
        rnms[SD6[[2]]] <- "DE"
        rownames(Simulated_Data_EE) <- rnms
        pe_mat <- rbind(pe_mat, Simulated_Data_EE)
        fcs <- c(fcs, rep(NA, nEE))
    }
    
    names(fcs) <- rownames(pe_mat)
    
    if (nDE + nDP + nDM + nDB + nEE + nEP == 0)
        stop("Please specify a nonzero number of 
            either DE, DP, DM, DB, EE, or EP genes.")
    
    SD <- list(Simulated_Data=pe_mat, FC=fcs)
    condition <- c(rep(1, numSamples), rep(2, numSamples))
    fcs <- data.frame(Category=names(fcs),
        FC=as.numeric(fcs))
    rownames(SD[[1]]) <- paste0(rownames(SD[[1]]), 1:nrow(SD[[1]]), sep="")
    rownames(fcs) <- rownames(SD[[1]])
    colnames(SD[[1]]) <- names(condition) <- paste0("Sample", 1:ncol(SD[[1]]), sep="")
    sce <- SingleCellExperiment(
        assays=list(counts=SD[[1]]), 
        colData=data.frame(condition),
        rowData=fcs)
    return(sce)
}


