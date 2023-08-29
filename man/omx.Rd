\name{omx}
\alias{omx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute full information maximum likelihood estimates of accuracy and precision for the
nonconstant bias measurement error model using 'OpenMx'. }
\description{
Compute full information maximum likelihood (FIML) estimates of accuracy (bias) and precision (imprecision) for the
nonconstant bias measurement error model.  'OpenMx' functions are used to construct and fit a one latent factor model. Likelihood-based confidence intervals and bootstrapped confidence interval are can be determined for all model parameters and relevant functions of the model parameters.
}
\usage{
omx(data, rvEst=rep(1,ncol(data)), mubarEst=mean(data[,1]), interval=0.95, 
  reps=500,bs.q=c(0.025,0.975), bs=TRUE)
}

\arguments{
  \item{data}{\eqn{n} (no. of items) \eqn{\times N} (no. of methods) \code{matrix} or \code{data.frame} containing the measurements. \eqn{N} must be \eqn{\ge 3}. Missing values are allowed.}
  \item{rvEst}{A \code{vector} of \eqn{N} residual variance starting values, one for each corresponding method.}
  \item{mubarEst}{A scalar starting value for estimating the the true mean value. }
  \item{interval}{Confidence level for likelihood-based confidence intervals. Should be a scalar value greater than 0 and less than 1.}
  \item{reps}{Number of bootstrap samples. Ignored if \code{bs=FALSE}.}
  \item{bs.q}{A \code{vector} of desired quantiles for bootstrapped samples. Default is \code{ci.q=c(0.025,0.975)}.}
  \item{bs}{A boolean indicating whether bootstapped samples are to be generated. Default is \code{TRUE}.}

}
\details{
Measurement Error Model:
  
\deqn{x_{ik}=\alpha_i+\beta_i\mu_k+\epsilon_{ik}}{x[i,k] = alpha[i] + beta[i]*mu[k] + epsilon[i,k]}
 
where \eqn{x_{ik}} is the measurement by the ith of \eqn{N} methods for the kth of \eqn{n} items, \eqn{i = 1} to \eqn{N\ge 3}, \eqn{k = 1} to \eqn{n}, \eqn{\mu_k} is the true value for the kth item, \eqn{\epsilon_{ik}} is the
normally distributed random error with variance \eqn{\sigma_i^2} for the ith method and the kth item, and
\eqn{\alpha_i} and \eqn{\beta_i} are the accuracy parameters for the ith method.  The beta for the first column of \code{data}) is set to one.  The corresponding alpha is set to 0.  These constraints or similar are required for model identification.

The imprecision for the ith method is \eqn{\sigma_i}. If all alphas are zeroes and all betas are ones, there is
no bias.  If all betas equal 1, then there is a constant bias.  If some of the betas differ from one there is a nonconstant bias.  Note that the individual betas are not unique - only ratios of the betas are unique. If you divide all the betas by \eqn{\beta_i}, then the betas represent the scale bias of the other devices/methods relative to device/method \eqn{i}.  Also, when the betas differ from one, the sigmas are not directly comparable because the measurement scales (size of the units) differ.  To make the sigmas comparable, divide them by their corresponding beta.

Technically, the alphas and betas describe the measurements in terms of the unknown true values (i.e., the unknown true values can be thought of as a latent variable).  The "true values" are ALWAYS unknown (unless you have a real, highly accurate reference method/device). The real goal is to calibrate one device/method in terms of another.  This is easily accomplished because each measurement is a linear function of the same unknown true values.  For methods 1 and 2, the calibration curve is given by:

\deqn{E[x_{1k}]=\left(\alpha_1-\alpha_2\beta_1/\beta_2\right)+\left(\beta_1/\beta_2\right)E[x_{2k}]}{E[x[1,k]] = (alpha[1] - alpha[2]*beta[1]/beta[2]) + (beta[1]/beta[2])*E[x[2,k]]}

or equivalently

\deqn{E[x_{2k}]=\left(\alpha_2-\alpha_1\beta_2/\beta_1\right)+\left(\beta_2/\beta_1\right)E[x_{1k}]}{E[x[2,k]] = (alpha[2] - alpha[1]*beta[2]/beta[1]) + (beta[2]/beta[1])*E[x[1,k]]}.

Use \code{\link{cplot}}, with the alpha.beta.sigma argument specified, to display this calibration curve, calibration equation, and the corresponding scale-bias adjusted imprecision standard deviations.

Note that likelihood confidence intervals and bootstrapped confidence intervals can be returned.  Wald-type intervals based on the standard errors are alos available by using the \code{confint} function on the returned \code{fit} object. See examples.

  
}
\value{
  \item{fit}{'OpenMx' fit object containing the results (FIML parameter estimates, etc) }
  \item{ci}{Likelihood-based confidence intervals for all parameters and certain useful functions of parameters.}
  \item{boot}{Object created by 'mxBootstrap' 'OpenMx' function. Not returned if \code{bs=FALSE}.}
  \item{q.boot}{\code{data.frame} containing the standard error and quantile estimates based on bootstrapped samples. Not returned if \code{bs=FALSE}.}
  \item{abs}{A \eqn{3 \times N} \code{matrix} of the estimated alphas, betas, and the raw imprecision standard deviations for each of the \eqn{N} methods. The results can be passed to the \code{merror} \code{cplot} to produce a calibration plot.}
  \item{bsReps}{Number of bootstrapped samples. Default is 500. Not returned if \code{bs=FALSE}.}
  \item{model}{The 'OpenMx' one-factor model. }
}

\note{

The following names are used to describe the estimates:

1) a1, a2, a3 and so forth denote the alphas (intercept).

2) b1, b2, b3 and so forth denote the betas (scale or slope).

3) ve1, ve2, ve3 and so forth denote the raw (uncorrected for scale) residual random error variances (imprecision variances).

4) se1 denotes the imprecison standard deviation for the reference method.

5) base2, base3 and so forth denote the scale bias-adjusted imprecision standard deviations on the scale of the reference method.

6) mubar is the estimated mean of the true values on the scale of the reference method.

7) sigma2 is the estimated variance of the true values on the scale of the reference method.

8) sigma is the estimated standard deviation of the true values on the scale of the reference method.
}

\author{ Richard A. Bilonick }

\seealso{ \code{\link{ncb.od}}, \code{\link{alpha.beta.sigma}} }
\examples{

\dontrun{
library(OpenMx)
library(merror)

data(pm2.5)

pm <- pm2.5

# OpenMx does not like periods in data column names
names(pm) <- c('ms_conc_1','ws_conc_1','ms_conc_2','ws_conc_2','frm')

# Fit model with FRM sampler as reference
omxfit <- omx(data=pm[,c(5,1:4)],bs.q=c(0.025,0.5,0.975),reps=100)

# Look at results
summary(omxfit$fit)$parameters[,c(1,5,6)] # Parameter estimates and standard errors
round(omxfit$ci[,1:3],3) # Likelihood-based intervals

# Estimated standard errors and quantiles based on bootstrapped samples
round(omxfit$q.boot,3)

# Wald-type intervals
#   - note variances not standard deviations and different ordering
confint(omxfit$fit) 

omxfit$abs # Use with cplot

# Make a calibration plot using the results from omx instead of the default ncb.od
cplot(pm[,c(5,1:4)],1,2,alpha.beta.sigma=omxfit$abs)
}
}

\keyword{ model }

