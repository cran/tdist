tdist <- function(funx, dff, lambda = rep(1, length(dff)), funtype = 1, pts = 14)
{
    iserr <- 1;
    xfun <- numeric(0);
    yfun <- numeric(0);
    # Choose correct function type
    if(!((funtype == 0) | (funtype == 1) | (funtype == 2) | (funtype == 3) | (funtype == 4)))   
    {
        print('!!! Choose funtype 1 for CDF, 2 for PDF, 3 for INV, or 4 for CHF'); break
    }
    # Dimensions of lambda and df should be equal
    if(length(lambda) == 1)
    {
        # If lambda is scalar resize it to size of df
        lambda <- rep(lambda, length(dff))
    }
    if(length(lambda) == 0) lambda <- rep(1, length(dff))            
    if(length(lambda) != length(dff))
    {
        print('!!! lambda should have the same size as df'); break
    }
    # Exclude random variables with zero coefficients in lambda
    if(any(lambda))
    {
        kk <- index(lambda)[lambda != 0]
        lambda <- lambda[kk]
        dff <- dff[kk]
    }
    # Exclude random variables with zero or negative coefficients in df  
    kk <- index(dff)[dff > 0]
    if(length(kk) !=length(dff))
    {
        iserr <- iserr*0
        print('!!! Excluded random variables with zero or negative dff !')
        lambda <- lambda[kk]
        dff <- dff[kk]
    }
    # Set yfun=NaN if all coefficients in df and/or in lambda are zeros
    if((!all(lambda)) | (!all(dff)))
    {
        yfun <- rep(NaN, length(funx))
        xfun <- funx
        iserr <- iserr*0
        print('!!! yfun=NaN if all coefficients in df and/or in lambda are zeros !')
        return(list('yfun' <- yfun, 'xfun' <- xfun, 'iserr' <- iserr)); break
    }
    # Variables with df > 100 set to be standard normal random variables
    kk <- index(dff)[dff == Inf]
    if(length(kk) != 0) dff[kk] <- 0
    kk <- index(dff)[dff > 100]
    if(length(kk) != 0)
    {
        iserr <- iserr*0
        print('!!! Variables with df > 100 were set to be standard normal random variables !')
        dff[kk] <- 0
    }
    # Find norm of lambda and estimate the approximate 95% range
    lambda <- abs(lambda)
    quant <- c(9.1, 3.6, 2.8, 2.5, 2.3, 2.2, 2.1, 2.1, 2.1, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0)
    dfind <- minv(ceiling(dff), 15)
    kk <- index(dfind)[dfind == 0]
    dfind[kk] <- 15
    qq <- quant[dfind]
    nr <- sqrt(sp(lambda, lambda))
    xupp <- sp(lambda, qq) / nr 
    # Generate an error message if large elements in funx
    if(length(funx) != 0)
    {
        xmax <- max(funx / nr)
        pcrit <- .005 * (((xupp - 2) / 7.1)^5)
        if ((xmax > 70) | ((funtype == 3) & ((1 - max(funx, 1 - min(funx))) < pcrit)))
        {
            iserr <- iserr*0
            print('!!! Large elements in funx ! The result could be inappropriate !')
        }
    }
    # Calculate GaussianQuadratureWeights
    resgrule <- grule(pts)
    bp <- resgrule$bp
    wf <- resgrule$wf
    xmax <- 15
    NN <- xmax * 10
    per <- pi / xmax
    limits <- c(0, 3 * per / 40, per, NN * per)
    ms <- c(3, 3, NN-1)
    resgweights <- gweights(limits, ms, bp, wf)
    tt <- resgweights$tt
    weights <- resgweights$w
    # Evaluate the characteristic function at funx
    if(funtype == 4)
    {
        if(length(funx) == 0) funx <- tt
        restchfvw <- tchfvw(funx, dff, lambda)
        yfun <- restchfvw$chf
        xfun <- restchfvw$tt
    }
    if(funtype != 4)
    {
        chf <- tchfvw(tt, dff, lambda / nr)$chf
        wf <- weights * chf
    }
    # Evaluate the required function
    if ((funtype == 0) | (funtype == 1) | (funtype == 2))
    {
        symetry <- 0
        if(length(funx) == 0)
        {
            funx <- nr * seq(0, xupp, , 51)
            symetry <- 1
        }
        yfun <- tcdfpdf(funx, tt, wf, funtype, nr)
        xfun <- funx
        if((symetry == 1) & (funtype == 0))
        {
            xsize <- length(funx)
            xfun <- c(-funx[seq(xsize, 2, -1)], funx)
            cdf <- yfun[, 1]
            pdf <- yfun[, 2]
            yfun <- matrix(0,length(xfun),2)
            yfun[, 1] <- c(1 - cdf[seq(xsize, 2, -1)], cdf)
            yfun[, 2] <- c(pdf[seq(xsize, 2, -1)], pdf)
        }
        if((symetry == 1) & (funtype == 1))
        {
            xsize <- length(funx)
            xfun <- c(-funx[seq(xsize, 2, -1)], funx)
            cdf <- yfun
            yfun <- c(1 - cdf[seq(xsize, 2, -1)], cdf)
        }
        if((symetry == 1) & (funtype == 2))
        {
            xsize <- length(funx)
            xfun <- c(-funx[seq(xsize, 2, -1)], funx)
            pdf <- yfun
            yfun <- c(pdf[seq(xsize, 2, -1)], pdf)
        }
    }
    if(funtype == 3)
    {
        restinvvw <- tinvvw(funx, tt, wf, nr)
        yfun <- restinvvw$yfun
        xfun <- restinvvw$xfun
        err <- restinvvw$err
        iserr <- iserr * err
    }
    # Change the iserr to be logical output
    if(iserr == 0) iserr <- 1
    if(iserr != 0) iserr <- 0
    # Output
    list('yfun' = yfun, 'xfun' = xfun, 'iserr' = iserr)
}
