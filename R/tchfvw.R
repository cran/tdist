tchfvw <- function(tt, nu, l)
{
    chf <- rep(1, length(tt))
    kk <- index(tt)
    kkk <- index(nu)[nu>0]
    if(length(kkk) != 0)
        for(k in 1:length(kkk))
        {
            a <- nu[kkk[k]] / 2
            b1 <- abs(l[kkk[k]] * tt[kk] * sqrt(nu[kkk[k]]))
            b2 <- (abs(l[kkk[k]] * tt[kk]) * sqrt(nu[kkk[k]]))^(nu[kkk[k]] / 2)
            b3 <- gamma(nu[kkk[k]] / 2) * (2^(nu[kkk[k]] / 2 - 1))
            chf[kk] <- chf[kk] * besselK(b1, a) * b2 / b3
        }
    kkk <- index(nu)[nu==0]
    if(length(kkk) != 0)
        for(k in 1:length(kkk))
            chf[kk] <- chf[kk] * exp(-1 * sqr(l[kkk[k]] * tt[kk]) / 2) 
    list('chf' = chf, 'tt' = tt)
}
