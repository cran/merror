\name{precision.mle.ncb.od}
\alias{precision.mle.ncb.od}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Computes iterative approximation to mle precision estimates for nonconstant bias model using original data. }
\description{
  This is an internal function that computes iterative approximation to mle precision estimates for nonconstant bias model using original data.
}
\usage{
precision.mle.ncb.od(x, M = 20, beta.bars = beta.bar(x), jaech.errors = FALSE)
}

\arguments{
  \item{x}{ A matrix or numeric data.frame consisting of an n (no. of items) by N (no. of methods) matrix of measuremnts. N must be greater than 3 and n > N. }
  \item{M}{ Maximum no. of iterations for convergence. }
  \item{beta.bars}{ Estimates or hypothesized values for the betas. }
  \item{jaech.errors}{ TRUE replicates the minor error in Jaech's Fortran code to allow comparison with his
  examples. }
}
\details{

Provides iterative approximation to MLE precision estimates for NonConstant Bias model using Original Data. See Jaech, p. 185-186.
 
}
\value{

  \item{sigma2 }{Estimated squared imprecisions (variances) for methods.}
  \item{sigma.mu2 }{Estimated process variance.}
}
\references{ Jaech, J. L. (1985) \emph{Statistical Analysis of Measurement Errors}. New York: Wiley. }
\author{ Richard A. Bilonick rab@nauticom.net}

\seealso{ \code{\link{precision.grubbs.ncb.od}},\code{\link{precision.grubbs.cb.pd}} }

\keyword{ htest }
