\name{TmixDist}
\alias{tdist}
\title{Distribution of a linear combination of independent Student's
t-variables}
\description{Computes the distribution of a linear combination of independent
Student's t random variables (with small degrees of freedom, df <= 100) and/or
standard Normal Z random variables.}
\usage{tdist(funx, dff, lambda = rep(1, length(dff)), funtype = 1, pts = 14)}
\arguments{
\item{funx}{vector of function (funtype) input values. If
\code{funx = numeric(0)} then \code{xfun} is generated automatically.}
\item{dff}{vector of degrees of freedom of independent Student's t random
variables, use \code{Inf} as df's for the standard Normal random variables.}
\item{lambda}{vector of coefficients of the linear combination.}
\item{funtype}{default value is 1 (calculates the cdf). The following funtypes
are legible: \cr
0: tdist calculates cdf and pdf at once, \code{yfun = cbind(cdf, pdf)}. \cr
1: tdist calculates the cumulative distribution function, cdf at funx,
\code{yfun = cdf}. \cr
2: tdist calculates the probability density function, pdf at funx,
\code{yfun = pdf}. \cr
3: tdist calculates the quantile function, qf at funx, \code{yfun = qf}. \cr
4: tdist calculates the characteristic function, chf at funx,
\code{yfun = chf}.}
\item{pts}{number of pts for Gaussian Quadrature. \cr
By default \code{pts = 14}. For many practical purposes, fast and reasonably
precise results are for choice of pts as small as 3.}
}
\value{A list of
\item{yfun}{vector with calculated function values, the result depends on
funtype. If \code{funtype = 0}, yfun has two columns (cdf and pdf).}
\item{xfun}{vector of function input values. Typically \code{xfun = funx}.
If \code{funx = numeric(0)}, xfun is generated automatically.}
\item{iserr}{Error message. If \code{iserr = 1}, some problem has occured during
calculation, see the warning message. If \code{iserr = 0}, corret calculations.}
}
\references{
Witkovsky, V. (2001), On the exact computation of the density and of the
quantiles of linear combinations of t and F random variables. \emph{Journal of
Statistical Planning and Inference}, 94, 1-13.

Witkovsky, V. (2004), Matlab algorithm TDIST: The distribution of a linear
combination of Student's T random variables. \emph{COMPSTAT 2004, 16th Symposium
of IASC PRAGUE, August 23-27}, Physica-Verlag/Springer 2004, 1995-2002.
}
\author{Viktor Witkovsky \email{witkovsky@savba.sk},
\url{http://aiolos.um.savba.sk/~viktor/}. Rewritten from Matlab algorithm to R
by Alexander Savin \email{savin@savba.sk}.}
\seealso{\code{\link{tcdfpdf}}, \code{\link{tchfvw}}, \code{\link{tinvvw}}.}
\examples{
# Plot the cdf of the random variable T = t_1 + 2*t_2 + 3*t_3 + 4*t_4 + Z,
# where Z is a random variable with standard normal distribution
# and t_1, t_2, t_3 and t_4 are random variables with Student's t distribution
# with 1, 2, 3 and 4 degrees of freedom. The random variables are assumed
# to be stochastically independent.

funx = numeric(0)
dff = c(1, 2, 3, 4, Inf)
lambda = c(1, 2, 3, 4, 1)
res = tdist(funx, dff, lambda, 1)
plot(res$xfun, res$yfun, type = 'l')

###
# Plot the pdf of the random variable T = ( Z + t_1 + t_10 )/3,
# where Z is a random variable with standard normal distribution
# and t_1 and t_10 are random variables with Student's t distribution
# with 1 and 10 degrees of freedom. The random variables are assumed
# to be stochastically independent.

funtype = 2
funx = numeric(0)
dff = c(Inf, 1, 10)
lambda = c(1, 1, 1) / 3
pts = 6
res = tdist(funx, dff, lambda, funtype, pts)
plot(res$xfun, res$yfun, type = 'l')

###
# Calculate the quantiles (for given probabilities 0.9, 0.95, 0.99)
# of the random variable  T = ( t_1 + Z )/2, where Z is a random variable
# with standard normal distribution and t_1 is a random variable
# with Student's t distribution with 1 degree of freedom. The random variables
# are assumed to be stochastically independent.

prob = c(0.9, 0.95, 0.99)
quantiles = tdist(prob, c(1, Inf), 1/2, 3)$yfun
cbind(prob, quantiles)

###
# Plot the characteristic function of the random variable T = t_1 + t_2 + t_3,
# where t_1, t_2, and t_3 are random variables with Student's t distribution
# with 1, 2, and 3 degrees of freedom. The random variables are assumed to be
# stochastically independent.

chftt = tdist(numeric(0), c(1, 2, 3), c(1, 1, 1), 4)
chf = chftt[[1]]
tt = chftt[[2]]
plot(tt, chf, type = 'l')
}
\keyword{distribution}
