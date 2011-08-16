\name{process.sd}
\alias{process.sd}
\title{Compute process standard deviation}
\description{This function computes the process standard
deviation and is used internally by the function
precision.grubbs.ncb.od.}
\usage{
process.sd(x)
}
\arguments{
  \item{x}{ A matrix or numeric data.frame consisting of an n (no. of items) by N
  (no. of methods) matrix of measuremnts. }
}
\value{
A scalar containing the method of moments estimate of the process standard deviation.
}
\details{
The process standard deviation is the standard deviation of the
true values uncontaminated by measurement error. See Jaech, p. 185.
}
\references{ Jaech, J. L. (1985) \emph{Statistical Analysis of Measurement Errors}. New York: Wiley. }
\author{ Richard A. Bilonick }

\seealso{ \code{\link{precision.grubbs.ncb.od}} }
\examples{

data(pm2.5)
process.sd(pm2.5) # estimate of the sd of the "true values using the method of moments")

}

\keyword{ htest }
