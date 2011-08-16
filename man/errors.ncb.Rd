\name{errors.ncb}
\alias{errors.ncb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Extracts the estimated measurement errors assuming there is a
nonconstant bias and using the original data values. }
\description{

Extracts the estimated measurement errors assuming there is a
nonconstant bias and using the original data values.  
}
\usage{
errors.ncb(x)
}

\arguments{
  \item{x}{  A matrix or numeric data.frame consisting of an n (no. of items) by N
  (no. of methods) matrix of measuremnts. N must be >= 3 and n > N. }
}
\details{
Errors should have a zero mean and should be Normally distributed with constant variance for a given method.
}
\value{  
\item{errors.ncb }{n x N matrix of estimated errors.}
}
\references{ Jaech, J. L. (1985) \emph{Statistical Analysis of Measurement Errors}. New York: Wiley. }
\author{ Richard A. Bilonick }

\seealso{ \code{\link{cb.pd}}, \code{\link{ncb.od}},\code{\link{lrt}} }
\examples{

data(pm2.5)
errors.ncb(pm2.5)

}

\keyword{ htest }
