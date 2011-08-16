\name{lrt}
\alias{lrt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Likelihood ratio test for all betas equalling one. }
\description{
  Likelihood ratio test statistic - H0: all betas = one.
}
\usage{
lrt(x, M = 40)
}

\arguments{
  \item{x}{ n (no. of items) x N (no. of methods) matrix or data.frame containing the measurements. N must be greater than 3 and n > N.}
  \item{M}{ Maximum no. of iterations for convergence. }
}
\details{
See Jaech, pp. 204-205.
}
\value{
  \item{n.items }{No.of items.}
  \item{N.methods }{No. of methods.'}
  \item{beta.bars }{N vector of estimated betas.}
  \item{lambda }{Chi-square test statistic.}
  \item{df }{Degrees of freedom for the test (N-1).'}
  \item{p.value }{Empirical significance level for the observed test statistic.'}
}
\references{ Jaech, J. L. (1985) \emph{Statistical Analysis of Measurement Errors}. New York: Wiley. }
\author{ Richard A. Bilonick }

\seealso{ \code{\link{ncb.od}},\code{\link{cb.pd}},\code{\link{pm2.5}} }
\examples{

data(pm2.5)

lrt(pm2.5) # compare all 5 samplers (4 personal and 1 frm)

lrt(pm2.5[,1:4]) # compare only the personal samplers

stem(lrt(pm2.5)$beta.bars) # examine the estimated betas

}
\keyword{ htest }
