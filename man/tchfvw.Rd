\name{TmixDistIF}
\alias{tchfvw}
\title{The internal function for coputation of the characteristic function of a
linear combination of independent Student's t-variables}
\description{Computes characteristic function of a linear combination of
independent Student's t-variables with small degrees of freedom (df <= 100).}
\usage{
tchfvw(tt, nu, l)
}
\arguments{
\item{tt}{vector of real number in wich characteristic function is evaluated.}
\item{nu}{vector of degrees of freedom}
\item{l}{vector of coefficients}
}
\value{
A list of
\item{chf}{vector of values of characteristic function.}
\item{tt}{vector of real number in wich characteristic function is evaluated.}
}
\author{Viktor Witkovsky \email{witkovsky@savba.sk}. Rewritten from Matlab
algorithm to R by Alexander Savin \email{savin@savba.sk}.}
\seealso{\code{\link{tcdfpdf}}, \code{\link{tdist}}, \code{\link{tinvvw}}.}
\keyword{distribution}
