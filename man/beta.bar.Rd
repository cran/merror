\name{beta.bar}
\alias{beta.bar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Compute the estimates of betas. }
\description{
  This function is used internally to compute the estimates of betas.
}
\usage{
beta.bar(x)
}

\arguments{
  \item{x}{ A matrix or numeric data.frame consisting of an n (no. of items) by N
  (no. of methods) matrix of measuremnts. }
}
\details{
See Jaech, p. 184.
}
\value{
 A vector of length N (no. of methods) containing the estimates of beta.
}
\references{ Jaech, J. L. (1985) \emph{Statistical Analysis of Measurement Errors}. New York: Wiley. }
\author{ Richard A. Bilonick }

\seealso{ \code{\link{cb.pd}}, \code{\link{ncb.od}},\code{\link{lrt}} }
\examples{

data(pm2.5)
beta.bar(pm2.5) # estimate betas (accuracy parameter)

}

\keyword{ htest }
