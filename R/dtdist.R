dtdist <- function(x, dff, lambda = rep(1, length(dff)), pts = 14, log.d = FALSE)
{
    res <- tdist(x, dff, lambda, 2, pts)$yfun
    if(log.d)
        res <- log(res)
    res
}
