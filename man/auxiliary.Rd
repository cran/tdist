\name{Auxiliary}
\alias{index}
\alias{maxv}
\alias{minv}
\alias{sp}
\alias{sqr}
\title{Auxiliary function}
\description{
Auxiliary functions:

index --- generates the vector of indices

maxv --- vectorized max

minv --- vectorized min

sp --- scalar product

sqr --- power race 2 (squared vector)
}
\usage{
index(x)

maxv(x, m)

minv(x, m)

sp(x, y)

sqr(x)
}
\arguments{
\item{x, y}{vectors.}
\item{m}{a value (scalar).}
}
\value{
\item{index}{vector of indices.}
\item{maxv}{vectorized max.}
\item{minv}{vectorized min.}
\item{sp}{scalar product of two vectors.}
\item{sqr}{power race 2 (squared vector).}
}
\author{Alexander Savin \email{savin@savba.sk}}
\keyword{arith}
