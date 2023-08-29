\name{cplot}
\alias{cplot}
\title{ Scatter plot of observations for a pair of devices with calibration curve. }
\description{
 Creates a scatter plot for any pair of observations in the data.frame and 
 superimposes the calibration curve.
}
\usage{
cplot(df, i, j, leg.loc="topleft", regress=FALSE, lw=1, t.size=1, alpha.beta.sigma=NULL)
}
\arguments{
  \item{df}{ n (no. of items) x N (no. of methods) matrix or data.frame containing the measurements. N must be >= 3 and n > N. }
  \item{i}{ Select column i for device i. }
  \item{j}{ Select column j for device j not equal to i.}
  \item{leg.loc}{ Location of the legend. }
  \item{regress}{ If TRUE, add both naive regression lines (for comparison only). }
  \item{lw}{ Line widths.}
  \item{t.size}{ Text size.}
  \item{alpha.beta.sigma}{ By default, \code{cplot} computes the bias (alpha, beta) and
imprecision (sigma) estimates using ncb.od. You can override this by specifying a
3 x N matrix of values with alpha on the first row, beta on the second row, and
sigma on the third row, in the same order as the methods.}
  
}
\details{By default, \code{cplot} displays the corresponding calibration curve for
devices i and j based on the parameter estimates for alpha, beta, and sigma computed
using \code{ncb.od}. You can overide this calibration curve by providing argument
alpha.beta.sigma with different estimates. Both naive regression lines (device i 
regressed on device j, and device j regressed on device i) by setting "regress=TRUE". 
Note, however, that the calibration curve will fall somewhere between these two
regression lines, depending on the the ratio of the imprecision standard deviations
(sigmas). (This may not hold if there are missing measurement data values given that ordinary regression requires deleting any item with one or more missing values.)
}
\value{
Produces a scatter plot with the calibration curve and titles that includes the
calibration equation and the scale-bias adjusted imprecision standard deviations.
}
\references{ Jaech, J. L. (1985) \emph{Statistical Analysis of Measurement Errors}. New York: Wiley. }
\author{ Richard A. Bilonick }

\seealso{ \code{\link{merror.pairs} }}
\examples{

library(merror)
data(pm2.5)

# Make various calibration plots for pm2.5 measurements
par(mfrow=c(2,2))
cplot(pm2.5,2,1)
cplot(pm2.5,3,1)
cplot(pm2.5,4,1)
# Add the naive regression lines JUST for comparison
cplot(pm2.5,5,1,regress=TRUE,t.size=0.9)

# This is redundant but illustrates using the
# argument alpha.beta.sigma
a <- ncb.od(pm2.5)$sigma.table$alpha.ncb[1:5]
b <- ncb.od(pm2.5)$sigma.table$beta[1:5]
s <- ncb.od(pm2.5)$sigma.table$sigma[1:5]

alpha.beta.sigma <- t(data.frame(a,b,s))

cplot(pm2.5,2,1,alpha.beta.sigma=alpha.beta.sigma)
cplot(pm2.5,2,1,alpha.beta.sigma=alpha.beta.sigma,regress=TRUE)
data(pm2.5)

\dontrun{
# Use omx function to specify the data for alpha.beta.sigma
pm <- pm2.5

# omx uses OpenMx which does not like periods in data column names
names(pm) <- c('ms_conc_1','ws_conc_1','ms_conc_2','ws_conc_2','frm')

# Fit one-factor measurement error model with FRM sampler as reference
omxfit <- omx(data=pm[,c(5,1:4)],bs.q=c(0.025,0.5,0.975),reps=100)

# Make a calibration plot using the results from omx instead of the default ncb.od
cplot(pm[,c(5,1:4)],1,2,alpha.beta.sigma=omxfit$abs)
}


}

\keyword{ htest }


