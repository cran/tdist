qtdist <- function(p, dff, lambda, pts =14, lower.tail = TRUE, log.p = FALSE)
{
    if(log.p)
        p <- exp(p)
    if(!lower.tail)
        p <- 1 - p
    res <- tdist(p, dff, lambda, 3, pts)$yfun
    res
}
