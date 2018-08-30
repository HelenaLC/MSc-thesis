tperm.fd <- function(x1fd,x2fd,nperm=200,q=0.05,argvals=NULL,plotres=TRUE,...) # first and second 
{                                                                          # groups of data,
    if( !is.fd(x1fd) | !is.fd(x2fd) ){                                     # number permuts
        stop("x1fd and x2fd must both be functional data objects")         # quantile
    }                                                                      # where to evaluate
                                                                           # do I plot
    rangeobs = x1fd$basis$range
    rangehat = x2fd$basis$range


    if( !prod(rangeobs == rangehat) ){
        stop("x1fd and x2fd do not have the same range.")
    }

    if(is.null(argvals)){
        argvals = seq(rangeobs[1],rangeobs[2],length.out=101)
    }

    q = 1-q

    x1mat = fda:::eval.fd(argvals,x1fd)
    x2mat = fda:::eval.fd(argvals,x2fd)

    n1 = ncol(x1mat)
    n2 = ncol(x2mat)

    Xmat = cbind(x1mat,x2mat)

    perms <- gtools::permutations(n1+n2, n1+n2, repeats.allowed=FALSE)
    g1 <- combn(n1+n2, n1)
    g2 <- apply(g1, 2, function(n) setdiff(seq_len(n1+n2), n))
    perms <- t(rbind(g1, g2))
    if (nperm > nrow(perms)) nperm <- nrow(perms)
    perms <- perms[sample(nperm), ]

    Tnullvals <- t(apply(perms, 1, function(idx) {
        tXmat = Xmat[, idx]
        
        tmean1 = rowMeans(tXmat[,1:n1])
        tmean2 = rowMeans(tXmat[,n1+(1:n2)])
        
        tvar1 = matrixStats::rowVars(tXmat[,1:n1])/n1
        tvar2 = matrixStats::rowVars(tXmat[,n1+(1:n2)])/n2
        
        denom <- sqrt(tvar1+tvar2)
        if (denom == 0) return(NA)
        abs(tmean1-tmean2)/denom
    }))
    Tnull <- matrixStats::colMaxs(Tnullvals, na.rm=TRUE)
    Tnull[is.infinite(Tnull)] <- NA

    mean1 = rowMeans(Xmat[,1:n1])
    mean2 = rowMeans(Xmat[,n1+(1:n2)])

    var1 = matrixStats::rowVars(Xmat[,1:n1])/n1
    var2 = matrixStats::rowVars(Xmat[,n1+(1:n2)])/n2

    Tvals = abs(mean1-mean2)/sqrt(var1+var2)
    Tvals[is.infinite(Tvals)] <- NA
    Tobs = max(Tvals, na.rm=TRUE)

    test = Tobs < Tnull
    test[is.na(test)] <- FALSE
    pval = mean(test, na.rm=TRUE)
    qval = quantile(Tnull, na.rm=TRUE)

    pvals.pts = rowMeans(Tvals < Tnullvals, na.rm=TRUE)
    qvals.pts = matrixStats::rowQuantiles(Tnullvals, probs=q)

    if(plotres){

	  if( is.null(names(x1fd$fdnames)) | is.null(names(x2fd$fdnames)) ){
		xlab='argvals'
	  }	
	  else if( prod(names(x1fd$fdnames)[1] == names(x2fd$fdnames)[1]) ){
		xlab = names(x1fd$fdnames)[1]
	  }
	  else{ xlab = 'argvals' }

        ylims = c( min(Tvals,qvals.pts),max(Tobs,qval))

        plot(argvals,Tvals,type='l',col=2,ylim=ylims,lwd=2,
		xlab=xlab,ylab='t-statistic',...)
        lines(argvals,qvals.pts,lty=3,col=4,lwd=2)
        abline(h=qval,lty=2,col=4,lwd=2)
	
        legendstr = c('Observed Statistic',
			    paste('pointwise',1-q,'critical value'),
			    paste('maximum',1-q,'critical value'))

	  legend(argvals[1],ylims[2],legend=legendstr,col=c(2,4,4),
		lty=c(1,3,2),lwd=c(2,2,2))
    }


    return( list(pval=pval,qval=qval,Tobs=Tobs,Tnull=Tnull,
        Tvals=Tvals,Tnullvals=Tnullvals,qvals.pts=qvals.pts,
        pvals.pts=pvals.pts,argvals=argvals) )
}
