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
  \item{x}{ n (no. of items) x N (no. of methods) \code{matrix} or \code{data.frame} containing the measurements. N must be >= 3. Missing values are not allowed. }
  \item{beta}{ N vector of betas, either estimated by \code{beta.bar} function or hypothesized. }
  \item{M}{ Maximum number of iterations for convergence. }
  \item{conf.level}{ Chosen confidence level which must be greater than zero and less than one. }
}
\details{
Measurement Error Model:
  
x[i,k] = alpha[i] + beta[i]*mu[k] + epsilon[i,k]
 
where x[i,k] is the measurement by the ith method for the kth item, i = 1 to N,
k = 1 to n, mu[k] is the true value for the kth item, epsilon[i,k] is the
normally distributed random error with variance sigma[i] squared for the ith method and the kth item, and
alpha[i] and beta[i] are the accuracy parameters for the ith method.  The product of the betas is constrained to equal one (equivalently, the geometric average of the beta's is constrained to one).  When the betas are all equal to one, the average of the alphas equals zero (equivalently, the sum of the alphas is constrained to zero).

The imprecision for the ith method is sigma[i]. If all alphas are zeroes and all betas are ones, there is
no bias.  If all betas equal 1, then there is a constant bias. If some of the betas differ from one there is a nonconstant bias.  Note that the individual betas are not unique - only ratios of the betas are unique.  If you divide all the betas by beta_i, then the betas represent the scale bias of the other devices/methods relative to device/method i.  Also, when the betas differ from one, the sigmas are not directly comparable because the measurement scales (size of the units) differ.  To make the sigmas comparable, divide them by their corresponding beta.  This result is shown as bias.adj.sigma.

By using the original data values, the betas can be estimated and also the process variance, that is, the variance of the true values.

Technically, the alphas and betas describe the measurements in terms of the unknown true values (i.e., the unknown true values can be thought of as a latent variable).  The "true values" are ALWAYS unknown (unless you have a real, highly accurate reference method/device).  The real goal is to calibrate one device/method in terms of another.  This is easily accomplished because each measurement is a function of the same unknow true values.  By solving the measurement error model (in expectation) for mu and substituting, any two devices/methods i=1 and i=2 can be be related as:

E[x[1,k]] = alpha[1] - alpha[2]*beta[1]/beta[2] + beta[1]/beta[2]*E[x[2,k]]

or equivalently

E[x[2,k]] = alpha[2] - alpha[1]*beta[2]/beta[1] + beta[2]/beta[1]*E[x[1,k]].

Use \code{\link{cplot}} to display this calibration curve and the corresponding scale-bias adjusted imprecision standard deviations.

The \code{omx} function is to be preferred to \code{ncb.od}. \code{omx} can accomodate missing measurement data values and can provide both likelihood-based confidence intervals and bootstrapped intervals for all parameters and relevant functions of parameters.

  
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
  \item{bias.adj.sigma }{sigma adjusted for scale bias: sigma/beta.}
  \item{H }{N+1 symmetric H matrix (see p. 201, Jaech).}
  \item{errors.nb }{n x N matrix of estimated measurement errors for no bias model.}
  \item{errors.cb }{n x N matrix of estimated measurement errors for constant bias model.}
  \item{errors.ncb }{n x N matrix of estimated measurement errors for nonconstant bias model}
}
\references{ Jaech, J. L. (1985) \emph{Statistical Analysis of Measurement Errors}. New York: Wiley. }
\author{ Richard A. Bilonick }

\seealso{ \code{\link{cb.pd}}, \code{\link{lrt}} }
\examples{

library(merror)
data(pm2.5)
ncb.od(pm2.5)               # nonconstant bias model using original data values
ncb.od(pm2.5,beta=rep(1,5)) # constant bias model using original data values

}

\keyword{ model }
