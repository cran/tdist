tcdfpdf <- function(x, tt, wf, funtype = 0, nr)
{
    # Integration by Gauss-quadrature
    x <- x / nr
    xsize <- length(x)
    pdf <- 0; cdf <- 0
    if(funtype == 0)
    {
        wft= wf / tt
        for(i in 1:xsize)
        {
            pdf[i]=sp(wf, cos(x[i] * tt))
            cdf[i]=sp(wft, sin(x[i] * tt))
        }
        pdf <- maxv(pdf / pi, 0) / nr
        cdf <- minv(maxv(1 / 2 + cdf / pi, 0), 1)
        yfun <- cbind(cdf, pdf) 
    }
    if(funtype == 1)
    {
        wft <- wf / tt
        for (i in 1:xsize) cdf[i] <- sp(wft, sin(x[i] * tt))
        cdf <- minv(maxv(1 / 2 + cdf / pi, 0), 1)
        yfun <- cdf
    }
    if(funtype == 2)
    {
        for(i in 1:xsize) pdf[i]=sp(wf, cos(x[i] * tt))
        pdf <- maxv(pdf / pi, 0) / nr
    yfun <- pdf
    }
    yfun
}
