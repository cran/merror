\name{alpha.beta.sigma}
\alias{alpha.beta.sigma}
\title{Build an alpha-beta-sigma Matrix for Use with the cplot Function}
\description{
  Creates a \eqn{3 \times N} (no. of methods) \code{matrix} consisting of the estimated alphas, betas, and imprecision sigmas for use with the \code{cplot} function.
}
\usage{
alpha.beta.sigma(x)
}
\arguments{
  \item{x}{A \eqn{k \times 3} \code{data.frame} with parameter estimates in the second column where \eqn{k} is the number of methods \eqn{m \times 3}.  The estimates should be arranged with the estimated \eqn{m-1} betas first, followed by the m residual variances, the variance of the true values, the \eqn{m-1} alphas, the mean of the true values. The \code{omx} function returns the fitted model in \code{fit} from which parameter estimates can be retrieved.  See the examples below.}
}

\details{
 This is primarily a helper function used by the \code{omx} function.
}

\value{
A \eqn{3 \times N} \code{matrix} consisting of alphas on the first row, betas on the second row, followed by raw imprecision sigmas.
}

\seealso{
  \code{\link{cplot}}, \code{\link{omx}}.
}
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

# Extract the estimates
alpha.beta.sigma(summary(omxfit$fit)$parameters[,c(1,5,6)])

# Make a calibration plot
cplot(pm[,c(5,1:4)],1,2,alpha.beta.sigma=
  alpha.beta.sigma(summary(omxfit$fit)$parameters[,c(1,5,6)]))

# The easier way
cplot(pm[,c(5,1:4)],1,2,alpha.beta.sigma=omxfit$abs)
}
}

\keyword{ model }