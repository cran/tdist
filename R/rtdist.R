rtdist <- function(n, dff, lambda)
{
    # Dimensions of lambda and dff should be equal
    if(length(lambda) == 1)
    {
        # If lambda is scalar resize it to size of dff
        lambda <- rep(lambda, length(dff))
    }
    if(length(lambda) == 0)
        lambda <- rep(1, length(dff))            
    if(length(lambda) != length(dff))
    {
        print('!!! lambda should have the same size as dff')
        break
    }
    # Exclude random variables with zero coefficients in lambda
    if(any(lambda))
    {
        kk <- index(lambda)[lambda != 0]
        lambda <- lambda[kk]
        dff <- dff[kk]
    }
    # Exclude random variables with zero or negative coefficients in dff
    kk <- index(dff)[dff > 0]
    if(length(kk) != length(dff))
    {
        print('!!! Excluded random variables with zero or negative dff !')
        lambda <- lambda[kk]
        dff <- dff[kk]
    }
    # Set yfun = NaN if all coefficients in df and/or in lambda are zeros
    if((!all(lambda)) | (!all(dff)))
    {
        yfun <- rep(NaN, length(funx))
        xfun <- funx
        print('!!! yfun=NaN if all coefficients in dff and/or in lambda are zeros !')
        return(list('yfun' <- yfun, 'xfun' <- xfun, 'iserr' <- iserr))
        break
    }
    # Variables with dff == Inf set to be standard normal random variables
    kk <- index(dff)[dff == Inf]
    if(length(kk) != 0) dff[kk] <- 0
    kk <- index(dff)[dff > 100]
    # If length(n) > 1, the length is taken to be the number required
    if(length(n) > 1)
        n <- length(n)
    res <- rep(0, n)
    for(i in 1:length(dff))
    {
       if(dff[i] == 0)            
       res <- res + lambda[i] * rnorm(n)
       if(dff[i] != 0)            
       res <- res + lambda[i] * rt(n, dff[i])
    }
    res
}
