ptdist <- function(q, dff, lambda, pts =14, lower.tail = TRUE, log.p = FALSE)
{
    if(lower.tail)
        res <- tdist(q, dff, lambda, 1, pts)$yfun
    if(!lower.tail)
        res <- 1 - tdist(q, dff, lambda, 1, pts)$yfun
    if(log.p)
        res <- log(res)
    res
}
