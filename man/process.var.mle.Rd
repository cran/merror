\name{process.var.mle}
\alias{process.var.mle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Compute process variance. }
\description{
  This is an internal function to compute the process variance.
}
\usage{
process.var.mle(sigma2, s, beta.bars, N, n)
}

\arguments{
  \item{sigma2}{ Estimated imprecisions for each method in an N vector. }
  \item{s}{ Variance-covariance N x N matrix. }
  \item{beta.bars}{ Estimates or hypothesized values for the N betas. }
  \item{N}{ No. of methods. }
  \item{n}{ No. of items. }
}
\details{
See Jaech p. 186  equations 6.37 - 6.3.10.
}
\value{

Estimated process variance.

}
\references{ Jaech, J. L. (1985) \emph{Statistical Analysis of Measurement Errors}. New York: Wiley. }
\author{ Richard A. Bilonick rab@nauticom.net}

\keyword{ htest }

