\name{mle}
\alias{mle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Compute maximum likelihood estimates of precision. }
\description{
  This is an internal function that computes the maximum likelihood estimates of precision for the constant
  bias model using paired data.
}
\usage{
mle(v, r, ni)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{v}{ Variance-Covariance matrix for the n x N items by methods measurement data. }
  \item{r}{ Initial estimates of imprecision, usually Grubbs. }
  \item{ni}{ No. of items measured. }
}

\value{
An N vector containing the maximum likelihood estimates of precision.
}
\references{ Jaech, J. L. (1985) \emph{Statistical Analysis of Measurement Errors}. New York: Wiley. }
\author{ Richard A. Bilonick }

\keyword{ htest }
