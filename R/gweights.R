gweights <- function(limits, ms, bp, wf)
{
    nl <- length(limits)
    tt <- numeric(0); w <- numeric(0)
    for(i in 1:(nl-1))
    {
        lower <- limits[i]
        upper <- limits[i+1]
        mparts <- ms[i]
        d <- (upper - lower) / mparts
        nquad <- length(bp)
        om <- rep(1, mparts)
        on <- rep(1, nquad)
        shift <- 1:mparts - 1
        ti <- d * (bp%*%t(om) + on%*%t(shift)) + lower
        wi <- d * wf%*%t(om)
        tt <- c(tt, as.vector(ti))
        w <- c(w, as.vector(wi))
    }
    list('tt' = tt, 'w' = w)
}
