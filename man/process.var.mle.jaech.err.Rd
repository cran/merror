\name{process.var.mle.jaech.err}
\alias{process.var.mle.jaech.err}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Compute process variance but with minor error in Jaech Fortran code. }
\description{
This is an internal function to compute the process variance that replicates the minor
error in Jaech's Fortran code. This allows comparing merror estimates to those
shown in Jaech 1985.
}
\usage{
process.var.mle.jaech.err(sigma2, s, beta.bars, N, n)
}

\arguments{
  \item{sigma2}{ Estimated imprecisions for each method in an N vector. }
  \item{s}{ Variance-covariance N x N matrix. }
  \item{beta.bars}{ Estimates or hypothesized values for the N betas }
  \item{N}{ No. of methods. }
  \item{n}{ No. of items. }
}
\details{
See Jaech p. 186  equations 6.37 - 6.3.10.
Jaech p. 288 line 2330 has s[i,j] instead of s[j,j].
Jaech p. 288 line 2410 omits "- 1/d2".  

}
\value{

Estimated process variance but replicating minor error in Jaech's Fortran code.
}
\references{ Jaech, J. L. (1985) \emph{Statistical Analysis of Measurement Errors}. New York: Wiley. }
\author{ Richard A. Bilonick rab@nauticom.net }

\keyword{ htest }
