\name{precision.grubbs.ncb.od}
\alias{precision.grubbs.ncb.od}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{  Computes Grubbs' method of moments estimators of precision for the
nonconstant bias model using original data values. }
\description{
This is an internal function that computes Grubbs' method of moments estimators of precision for the
nonconstant bias model using original data values.
}
\usage{
precision.grubbs.ncb.od(x, beta.bar.x = beta.bar(x))
}

\arguments{
  \item{x}{ A matrix or numeric data.frame consisting of an n (no. of items) by N
  (no. of methods) matrix of measuremnts. N must be >= 3 and n > N. }
  \item{beta.bar.x}{ Either estimates of beta or hypothesized values (one for
  each method in an N vector). }
}
\details{
See Jaech, p. 184.
}
\value{
Grubbs' method of moments estimates of the squared imprecision (variances).
}
\references{ Jaech, J. L. (1985) \emph{Statistical Analysis of Measurement Errors}. New York: Wiley. }
\author{ Richard A. Bilonick }
 
\seealso{\code{\link{precision.grubbs.cb.pd}}, \code{\link{ncb.od}}, \code{\link{cb.pd}},\code{\link{lrt}} }

\keyword{ htest }
