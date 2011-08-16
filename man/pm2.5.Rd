\name{pm2.5}
\alias{pm2.5}
\docType{data}
\title{ PM 2.5 Concentrations from SCAMP Collocated Samplers}
\description{
  Five filter-based samplers for measuring PM 2.5 concentrations were collocated and
  provided 77 complete sets of concentrations. This data was collected by the
  Stuebenville Comprehensive Air Monitoring Program (SCAMP) to check the accuracy and
  precision of the instruments.
}
\usage{data(pm2.5)}
\format{
  A data frame with 77 sets of PM 2.5 concentrations (micrograms per cubic meter) from the following 5 samplers:
  \describe{
    \item{ms.conc.1}{- personal sampler 1 - filter MS}
    \item{ws.conc.1}{- personal sampler 1 - filter WS}
    \item{ms.conc.2}{- personal sampler 2 - filter MS}
    \item{ws.conc.2}{- personal sampler 2 - filter WS}
    \item{frm}{- Federal Reference Method sampler}
  }
}
\source{
 Stuebenville Comprehensive Air Monitoring Program (SCAMP) 
}
\examples{
data(pm2.5)
boxplot(pm2.5)
merror.pairs(pm2.5)

# estimates of accuracy and precision
#   for nonconstant bias model using 
#   original data values
ncb.od(pm2.5)
}
\keyword{ datasets }

