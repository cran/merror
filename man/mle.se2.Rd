\name{mle.se2}
\alias{mle.se2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Compute squared standard errors for imprecision estimates for the
constant bias model using paired data. }
\description{
  This is an internal function that computes squared standard errors for imprecision estimates of the
constant bias model using paired data.
}
\usage{
mle.se2(v, r, ni)
}

\arguments{
  \item{v}{ Variance-Covariance matrix for the n x N items by methods measurement data. }
  \item{r}{ Initial estimates of imprecision, usually Grubbs }
  \item{ni}{ No. of items measured }
}
\details{
Computes the squared standard errors for the squared precisions. Before calling this function, compute the MLE's
}
\value{
  An N+1 symmetric H matrix. See p. 201 of Jaech, 1985, eq. 6.4.2.
}
\references{ J. L. Jaech, Statistical Analysis of Measurement Errors, Wiley, New
York: 1985. }
\author{ Richard A. Bilonick rab@nauticom.net }

\keyword{ htest }
