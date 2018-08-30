simulateSet <- function(SCdat, index, numSamples = 100,
    nDE = 250, nDP = 250, nDM = 250, nDB = 250, nEE = 5000, nEP = 4000,
    sd.range = c(1, 3), modeFC = c(2, 3, 4), random.seed = 284,
    varInflation = NULL, condition = "condition",
    param = BiocParallel::bpparam()) {

    ref <- unique(colData(SCdat)[[condition]])[1]
    BiocParallel::register(BPPARAM = param)
    FC <- scDD:::findFC(SCdat, index, sd.range = sd.range, N = 6,
        overExpressionProb = 0.5, plot.FC = FALSE, condition)
    constantZero <- NULL
    generateZero <- "empirical"
    conds <- colData(SCdat)[[condition]]
    Dataset1 <- assays(SCdat[, conds == ref])$normcounts
    set.seed(random.seed)
    pe_mat <- NULL
    fcs <- NULL
    if (nDE > 0) {
        SD1 <- scDD:::singleCellSimu(Dataset1, Method = "DE", index,
            FC, modeFC, Validation = FALSE, numGenes = nDE, numDE = nDE,
            numSamples = numSamples, generateZero = generateZero,
            constantZero = constantZero, varInflation)
        Simulated_Data <- SD1[[1]]
        rnms <- rep("EE", nrow(Simulated_Data))
        rnms[SD1[[2]]] <- "DE"
        rownames(Simulated_Data) <- rnms
        pe_mat <- Simulated_Data
        fcs <- SD1[[3]]
    }
    if (nDP > 0) {
        SD2 <- scDD:::singleCellSimu(Dataset1, Method = "DP", index,
            FC, modeFC, DP = c(0.33, 0.66), Validation = FALSE,
            numGenes = nDP, numDE = nDP, numSamples = numSamples,
            generateZero = generateZero, constantZero = constantZero,
            varInflation = varInflation)
        Simulated_Data_DP <- SD2[[1]]
        rnms <- rep("EE", nrow(Simulated_Data_DP))
        rnms[SD2[[2]]] <- "DP"
        rownames(Simulated_Data_DP) <- rnms
        pe_mat <- rbind(pe_mat, Simulated_Data_DP)
        fcs <- c(fcs, SD2[[3]])
    }
    if (nDM > 0) {
        SD3 <- scDD:::singleCellSimu(Dataset1, Method = "DM", index,
            FC, modeFC, Validation = FALSE, numGenes = nDM, numDE = nDM,
            numSamples = numSamples, generateZero = generateZero,
            constantZero = constantZero, varInflation)
        Simulated_Data_DM <- SD3[[1]]
        rnms <- rep("EE", nrow(Simulated_Data_DM))
        rnms[SD3[[2]]] <- "DM"
        rownames(Simulated_Data_DM) <- rnms
        pe_mat <- rbind(pe_mat, Simulated_Data_DM)
        fcs <- c(fcs, SD3[[3]])
    }
    if (nDB > 0) {
        SD4 <- scDD:::singleCellSimu(Dataset1, Method = "DB", index,
            FC, modeFC, DP = c(0.5, 0.5), Validation = FALSE,
            numGenes = nDB, numDE = nDB, numSamples = numSamples,
            generateZero = generateZero, constantZero = constantZero,
            varInflation)
        Simulated_Data_DB <- SD4[[1]]
        rnms <- rep("EE", nrow(Simulated_Data_DB))
        rnms[SD4[[2]]] <- "DB"
        rownames(Simulated_Data_DB) <- rnms
        pe_mat <- rbind(pe_mat, Simulated_Data_DB)
        fcs <- c(fcs, SD4[[3]])
    }
    if (nEP > 0) {
        SD5 <- scDD:::singleCellSimu(Dataset1, Method = "DP", index,
            FC, modeFC, DP = c(0.5, 0.5), Validation = FALSE,
            numGenes = nEP, numDE = 0, numSamples = numSamples,
            generateZero = generateZero, constantZero = constantZero,
            varInflation)
        Simulated_Data_EP <- SD5[[1]]
        rnms <- rep("EP", nrow(Simulated_Data_EP))
        rnms[SD5[[2]]] <- "DP"
        rownames(Simulated_Data_EP) <- rnms
        pe_mat <- rbind(pe_mat, Simulated_Data_EP)
        fcs <- c(fcs, SD5[[3]])
    }
    if (nEE > 0) {
        SD6 <- scDD:::singleCellSimu(Dataset1, Method = "DE", index,
            FC, modeFC, Validation = FALSE, numGenes = nEE, numDE = 0,
            numSamples = numSamples, generateZero = generateZero,
            constantZero = constantZero, varInflation = varInflation)
        Simulated_Data_EE <- SD6[[1]]
        rnms <- rep("EE", nrow(Simulated_Data_EE))
        rnms[SD6[[2]]] <- "DE"
        rownames(Simulated_Data_EE) <- rnms
        pe_mat <- rbind(pe_mat, Simulated_Data_EE)
        fcs <- c(fcs, rep(NA, nEE))
    }
    names(fcs) <- rownames(pe_mat)
    SD <- list(Simulated_Data = pe_mat, FC = fcs)
    condition <- c(rep(1, numSamples), rep(2, numSamples))
    fcs <- data.frame(Category = names(fcs), FC = as.numeric(fcs))
    rownames(SD[[1]]) <- paste0(rownames(SD[[1]]), 1:nrow(SD[[1]]), sep = "")
    rownames(fcs) <- rownames(SD[[1]])
    colnames(SD[[1]]) <- names(condition) <-
        paste0("Sample", 1:ncol(SD[[1]]), sep = "")
    SingleCellExperiment(
        assays = list(counts = SD[[1]]),
        colData = data.frame(condition), rowData = fcs)
}
