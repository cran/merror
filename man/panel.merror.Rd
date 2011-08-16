\name{panel.merror}
\alias{panel.merror}
\title{Draw diagonal line (line of equality) on merror.pairs plots}
\description{This functionis used internally by the function
merror.pairs.}
\usage{
panel.merror(x,y, ...)
}
\arguments{
  \item{x}{ A vector of measurements for one device, of length n. }
  \item{y}{ A vector of measurements for another device, of length n. }
  \item{...}{Additional arguments.}
}
\value{
Draws the diagonal line that represents the "line of equality", i.e., the "no bias model".
}
\author{ Richard A. Bilonick }

\seealso{ \code{\link{merror.pairs}} }

\keyword{ htest }
