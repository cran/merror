\name{redshift}
\alias{redshift}
\docType{data}
\title{ Spectroscopic and Photometric Galaxy Redshift Measurements }
\description{
The redshift observations were taken from DEEP 2 Galaxy Redshift Survey.  
}
\usage{data(redshift)}
\format{

Redshift measurements are usually denoted by \eqn{z}. 

A data frame with one spectroscopic redshift measurement and six different photometric measurements (by researcher) for 1432 galaxies:
  \describe{
    \item{z_spec}{Spectroscopic redshift}
    \item{z_fink}{Photometric redshift - S. Finklestein}
    \item{z_font}{Photometric redshift - A. Fontana}
    \item{z_pfor}{Photometric redshift - J. Pforr}
    \item{z_salv}{Photometric redshift - M. Salvator}
    \item{z_wikl}{Photometric redshift - T. Wiklind}
    \item{z_wuyt}{Photometric redshift - S. Wuyts}
  }
}
\details{

Because the photometric methods depend on the same color information, a one-factor measurement error model incuding both the spectroscopic and photomentric measurements would not be a viable model because the photometric measurements would tend to be correlated. A two-factor model would be needed but would require at minimum replicated spectroscopic measurments.
}
\source{
 
<https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20140013340.pdf>

}
\references{
Newman, Jeffrey A., Michael C. Cooper, Marc Davis, S. M. Faber, Alison L. Coil, Puragra Guhathakurta, David C. Koo et al. "The DEEP2 Galaxy Redshift Survey: Design, observations, data reduction, and redshifts." The Astrophysical Journal Supplement Series 208, no. 1 (2013): 5.
}

\examples{

library(OpenMx)
library(merror)

data(redshift)
merror.pairs(redshift)

# estimates of accuracy and precision
#   parameters for a one-factor
#   measurement error model
head(redshift)
merror.pairs(redshift)

\dontrun{
red <- omx(redshift[,-1],reps=200) # Drop the spectroscopic measurements

summary(red$fit)
red$ci
red$q.boot

cplot(redshift[,-1],1,2,alpha.beta.sigma=red$abs)
}}
\keyword{ datasets }

