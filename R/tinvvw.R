tinvvw <- function(p, tt, wf, nr)
{
    err <- 1
    if(length(p) == 0) p <- seq(0, 1, , 51)
    qq <- rep(0, length(p))
    # For p=0 set the quantile to be -Inf
    k <- index(p)[p == 0]
    if(any(k))
    {
        tmp <- -1 * Inf
        qq[k] <- rep(tmp, length(k))
    }
    # For p=1 set the quantile to be Inf
    k <- index(p)[p==1]
    if(any(k))
    {
        tmp <- Inf
        qq[k] <- rep(tmp, length(k))
    }
    # Newton's method
    clim <- 100
    count <- 0
    k <- index(p)[p > 0 & p < 1]
    if(length(k) == 0) yfun <- qq
    if(length(k) != 0)
    {
        xfun <- p
        pk <- xfun[k];
        xk <- rep(0, length(pk))
        h <- rep(1, length(pk))
        crit <- 1e-12
        while(any(abs(h) > crit * abs(xk)) & (max(abs(h)) > crit) & (count < clim))
        {
            count <- count + 1;
            yy <- tcdfpdf(xk, tt, wf, 0, nr)
            h <- (yy[,1] - pk) / yy[,2]
            xk <- xk - h
        }
        qq[k] <- xk
        if(count == clim)
        {
            err=0
            print('!!! TINVVW did not converge')
        }
        yfun <- qq
    }
    list('yfun' = yfun, 'xfun' = xfun, 'err' = err)
}
