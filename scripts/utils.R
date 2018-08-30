# confusion matrix
get_status <- function(U, V) { 
    if (any(is.na(U))) { 
        # this will happen when a method fails to run
        return(rep(0, 4))
    } else {
        # V is truth!
        # TP: pair of points in same cluster in U, V
        # TN: pair of points in diff cluster in U, V
        # FP: pair of points in same cluster in U, diff in V
        # FN: pair of points in diff cluster in U, same in V
        # when there are NA values in the predictions: FN++
        bu <- outer(U, U, '==')
        bv <- outer(V, V, '==')
        TP <- length(intersect(which( bu), which( bv))) - length(U)
        TN <- length(intersect(which(!bu), which(!bv)))
        FP <- length(intersect(which( bu), which(!bv)))
        FN <- length(intersect(which(!bu), which( bv))) + sum(is.na(bu))
        return(c(TP,TN,FP,FN)/2)
    }
}

# clustering evaluation metrics
precision <- function(x) {
    c(TP,TN,FP,FN) %<-% x
    TP / (TP + FP)
}
recall <- function(x) {
    c(TP,TN,FP,FN) %<-% x
    TP / (TP + FN)
}
accuracy <- function(x) {
    c(TP,TN,FP,FN) %<-% x
    (TP + TN) / sum(x)
}
f1 <- function(x) {
    c(TP,TN,FP,FN) %<-% x
    p <- TP / (TP+FP)
    r <- TP / (TP+FN)
    2 / (1/p + 1/r)
}
