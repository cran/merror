\name{ncb.od}
\alias{ncb.od}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Compute accuracy estimates and maximum likelihood estimates of precision for the
nonconstant bias measurement error model using original data. }
\description{
 Compute accuracy estimates and maximum likelihood estimates of precision for the
nonconstant bias measurement error model using original data.
}
\usage{
ncb.od(x, beta = beta.bar(x), M = 40, conf.level = 0.95)
}
\arguments{
  \item{x}{ n (no. of items) x N (no. of methods) matrix or data.frame containing the measurements. N must be greater than 3 and n > N. }
  \item{beta}{ N vector of betas, either estimated by \code{beta.bar} function or hypothesized. }
  \item{M}{ Maximum no. of iterations for convergence. }
  \item{conf.level}{ Chosen confidence level. }
}
\details{
Measurement Error Model:
  
x[i,k] = alpha[i] + beta[i]*mu[k] + epsilon[i,k]
 
where x[i,k] is the measurement by the ith method for the kth item, i = 1 to N,
k = 1 to n, mu[k] is the true value for the kth item, epsilon[i,k] is the
Normally distributed random error with variance sigma[i] squared for the ith method and the kth item, and
alpha[i] and beta[i] are the accuracy parameters for the ith method.

The imprecision for the ith method is sigma[i]. If all alphas are zeroes and all betas are ones, there is
no bias. If all betas equal 1, then there is a constant bias. Otherwise there is a nonconstant bias.

By using the original data values, the betas can be estimated and also the process variance.
  
}
\value{
  \item{conf.level }{Confidence level used.}
  \item{sigma.table }{Table of accuracy and precision estimates and confidence
  intervals.}
  \item{n.items }{No. of items.}
  \item{N.methods }{No. of methods}
  \item{sigma2 }{N vector of variances that measure the method imprecision.}
  \item{alpha.cb }{N vector of estimated alphas for constant bias model.}
  \item{alpha.ncb }{N vector of estimated alphas for nonconstant bias model.}
  \item{beta }{N vector of estimated or hypothesized betas.}
  \item{df }{N vector of estimated degrees of freedom.}
  \item{lb }{N vector of lower bounds for confidence intervals.}
  \item{ub }{N vector of upper bounds for confidence intervals.}
  \item{H }{N+1 symmetric H matrix (see p. 201, Jaech).}
  \item{errors.nb }{n x N matrix of estimated measurement errors for no bias model.}
  \item{errors.cb }{n x N matrix of estimated measurement errors for constant bias model.}
  \item{errors.ncb }{n x N matrix of estimated measurement errors for nonconstant bias model}
}
\references{ Jaech, J. L. (1985) \emph{Statistical Analysis of Measurement Errors}. New York: Wiley. }
\author{ Richard A. Bilonick rab@nauticom.net }

\seealso{ \code{\link{cb.pd}},\code{\link{lrt}} }
\examples{

data(pm2.5)
ncb.od(pm2.5)               # nonconstant bias model using original data values
ncb.od(pm2.5,beta=rep(1,5)) # constant bias model using original data values

}

\keyword{ htest }
