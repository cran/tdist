\name{TmixDist}
\alias{ctdist}
\alias{dtdist}
\alias{ptdist}
\alias{qtdist}
\alias{rtdist}
\title{Distribution of a linear combination of independent Student's
t-variables}
\description{Characteristic function, density, distribution function, quantile
function and random generation for the linear combination of independent
Student's t random variables (with small degrees of freedom, df <= 100) and/or
standard Normal Z random variables.}
\usage{
ctdist(x, dff, lambda, pts =14)
dtdist(x, dff, lambda, pts =14, log.d = FALSE)
ptdist(q, dff, lambda, pts =14, lower.tail = TRUE, log.p = FALSE)
qtdist(p, dff, lambda, pts =14, lower.tail = TRUE, log.p = FALSE)
rtdist(n, dff, lambda)
}
\arguments{
\item{x,q}{vector of quantiles. If \code{x = numeric(0)}, \code{q = numeric(0)}
then \code{x}, \code{q} is generated automatically.}
\item{p}{vector of probabilities. If \code{p = numeric(0)} then \code{p} is
generated automatically.}
\item{n}{number of observations. If \code{length(n) > 1}, the length is taken to
be the number required.}
\item{dff}{vector of degrees of freedom of independent Student's t random
variables, use \code{Inf} as df's for the standard Normal random variables.}
\item{lambda}{vector of coefficients of the linear combination.}
\item{pts}{number of pts for Gaussian Quadrature. \cr
By default \code{pts = 14}. For many practical purposes, fast and reasonably
precise results are for choice of pts as small as 3.}
\item{log.d}{logical; if TRUE, \code{dtdist} gives the logarithm of density.}
\item{log.p}{logical; if TRUE, probabilities p are given as log(p).}
\item{lower.tail}{logical; if TRUE (default), probabilities are P[X <= x],
otherwise, P[X > x].}
}
%\details{}
\value{\code{ctdist} gives the characteristic function, \code{dtdist} gives the
density, \code{ptdist} gives the distribution function, \code{qtdist} gives the
quantile function, and \code{rdist} generates
random deviates.}
\references{Witkovsky, V. (2001), On the exact computation of the density and of
the quantiles of linear combinations of t and F random variables. \emph{Journal
of Statistical Planning and Inference}, 94, 1-13.

Witkovsky, V. (2004), Matlab algorithm TDIST: The distribution of a linear
combination of Student's T random variables. \emph{COMPSTAT 2004, 16th Symposium
of IASC PRAGUE, August 23-27}, Physica-Verlag/Springer 2004, 1995-2002.}
\author{Alexander Savin \email{savin@savba.sk} and Viktor Witkovsky
\email{witkovsky@savba.sk}, \url{http://aiolos.um.savba.sk/~viktor/}.}
\seealso{\code{\link{tdist}}}
\examples{
# Plot the characteristic function of the random variable T = t_1 + t_2 + t_3,
# where t_1, t_2, and t_3 are random variables with Student's t distribution
# with 1, 2, and 3 degrees of freedom. The random variables are assumed to be
# stochastically independent.

tt = c(seq(0,.1,,100), seq(.1,7,,200))
chf = ctdist(tt, c(1, 2, 3), c(1, 1, 1))
plot(tt, chf, type = 'l')

###
# Plot the pdf of the random variable T = ( Z + t_1 + t_10 )/3,
# where Z is a random variable with standard normal distribution
# and t_1 and t_10 are random variables with Student's t distribution
# with 1 and 10 degrees of freedom. The random variables are assumed
# to be stochastically independent.

dff = c(Inf, 1, 10)
lambda = c(1, 1, 1) / 3
x = seq(-4,4,,100)
pdf = dtdist(x, dff, lambda, 6)
plot(x, pdf, type = 'l')

###
# Plot the cdf of the random variable T = t_1 + 2*t_2 + 3*t_3 + 4*t_4 + Z,
# where Z is a random variable with standard normal distribution
# and t_1, t_2, t_3 and t_4 are random variables with Student's t distribution
# with 1, 2, 3 and 4 degrees of freedom. The random variables are assumed
# to be stochastically independent.

dff = c(1, 2, 3, 4, Inf)
lambda = c(1, 2, 3, 4, 1)
x = seq(-36,36,,200)
cdf = ptdist(x, dff, lambda)
plot(x, cdf, type = 'l')

###
# Calculate the quantiles (for given probabilities 0.9, 0.95, 0.99)
# of the random variable  T = ( t_1 + Z )/2, where Z is a random variable
# with standard normal distribution and t_1 is a random variable
# with Student's t distribution with 1 degree of freedom. The random variables
# are assumed to be stochastically independent.

prob = c(0.9, 0.95, 0.99)
quantiles = qtdist(prob, c(1, Inf), 1/2, 3)
cbind(prob, quantiles)

###
# Generate 1000 realization of the random variable  T = ( t_3 + t_5 + Z )/3,
# where Z is a random variable with standard normal distribution and t_3 and t_5
# are a random variable with Student's t distribution with 3 and 5 degrees
# of freedom. The random variables are assumed to be stochastically independent.

r = rtdist(1000, c(3, 5, Inf), 1/3)

# Plot the binned kernel density estimate of the probability density of the
# generated data by 'bdke' procedure and compare with true density (red line).

library(KernSmooth)
est <- bkde(r, bandwidth=0.25)
plot(est, ylim = c(-0, .6), type="l")
lines(est$x,dtdist(est$x, c(3, 5, Inf), 1/3), col = 'red')
}
\keyword{distribution}
