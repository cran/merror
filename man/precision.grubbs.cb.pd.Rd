\name{precision.grubbs.cb.pd}
\alias{precision.grubbs.cb.pd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Computes Grubbs' method of moments estimators of precision for the
constant bias model using paired differences. }
\description{
This is an internal function that computes Grubbs' method of moments estimators of precision for the
constant bias model using paired differences 
}
\usage{
precision.grubbs.cb.pd(x)
}
\arguments{
  \item{x}{ A matrix or numeric data.frame consisting of an n (no. of items) by N
  (no. of methods) matrix of measuremnts. N must be greater than 3 and n > N. }
}
\details{
See Jaech 1985, Chapters 3 & 4, p. 144 in particular.
}
\value{
Estimated squared imprecision estimates (variances).
}
\references{ Jaech, J. L. (1985) \emph{Statistical Analysis of Measurement Errors}. New York: Wiley.}
\author{ Richard A. Bilonick rab@nauticom.net}
 
\seealso{ \code{\link{precision.grubbs.ncb.od}}, \code{\link{ncb.od}}, \code{\link{cb.pd}},\code{\link{lrt}} }

\keyword{ htest }
