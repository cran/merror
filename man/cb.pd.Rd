\name{cb.pd}
\alias{cb.pd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Compute accuracy estimates and maximum likelihood estimates of precision for the
constant bias measurement error model using paired data. }
\description{
Compute accuracy estimates and maximum likelihood estimates of precision for the
constant bias measurement error model using paired data.
}
\usage{
cb.pd(x, conf.level = 0.95, M = 40)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ n (no. of items) x N (no. of methods) matrix or data.frame containing the measurements. N must be greater than 3 and n > N.}
  \item{conf.level}{ Chosen onfidence level. }
  \item{M}{ Maximum no.of iterations to reach convergence. }
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

ME (method of moments estimator) and MLE are the same for N=3 instruments
 except for a factor of (n-1)/n: MLE = (n-1)/n * ME

Using paired differences forces Constant Bias model (beta[1] = beta[2] = ... = beta[N]). Also,
the process variance CANNOT be estimated.
}
\value{
  \item{conf.level }{Confidence level used.}
  \item{sigma.table }{Table of accuracy and precision estimates and confidence
  intervals.}
  \item{n.items }{No. of items.}
  \item{N.methods }{No. of methods}
  \item{Grubbs.initial.sigma2 }{N vector of initial imprecision estimates using Grubbs' method}
  \item{sigma2 }{N vector of variances that measure the method imprecision.}
  \item{sigma2.se2 }{N vector of squared standard errors of the estimated imprecisions (variances).}
  \item{alpha.cb }{N vector of estimated alphas for constant bias model.}
  \item{alpha.ncb }{N vector of estimated alphas for nonconstant bias model}
  \item{beta }{N vector of hypothesized betas for the constant bias model - all ones.}
  \item{df }{N vector of estimated degrees of freedom.}
  \item{chisq.low }{N vector of chi-square values for the lower tail (used to compute the ci upper bound).}
  \item{chisq.low }{N vector of chi-square values for the upper tail (used to compute the ci lower bound).}
  \item{lb }{N vector of lower bounds for confidence intervals}
  \item{ub }{N vector of upper bounds for confidence intervals}
}
\references{ Jaech, J. L. (1985) \emph{Statistical Analysis of Measurement Errors}. New York: Wiley. }
\author{ Richard A. Bilonick rab@nauticom.net }

\seealso{ \code{\link{ncb.od}},
 \code{\link{lrt}}
 }
\examples{

data(pm2.5)
cb.pd(pm2.5)

}
\keyword{ htest }

