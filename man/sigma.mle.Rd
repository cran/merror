\name{sigma.mle}
\alias{sigma.mle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Computes the ith iteration for computing the
  squared imprecision estimates. }
\description{
  This is an internal function that computes the ith iteration for computing the
  squared imprecision estimates.
}
\usage{
sigma.mle(i, s, sigma2, sigma.mu2, beta.bars, N, n)
}

\arguments{
  \item{i}{ Iteration i. }
  \item{s}{ Variance-covariance N x N matrix. }
  \item{sigma2}{ Estimated imprecisions for each method in an N vector }
  \item{sigma.mu2}{ Estimated process varinace. }
  \item{beta.bars}{ Estimates or hypothesized values for the N betas. }
  \item{N}{ No. of methods. }
  \item{n}{ No. of items. }
}
\details{
See Jaech p. 185-186 equations 6.3.1 - 6.3.6.
}
\value{
Estimated squared imprecisions (variances) for the ith iteration.
}
\references{ Jaech, J. L. (1985) \emph{Statistical Analysis of Measurement Errors}. New York: Wiley. }
\author{ Richard A. Bilonick rab@nauticom.net }

\keyword{ htest }
