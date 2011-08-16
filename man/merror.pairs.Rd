\name{merror.pairs}
\alias{merror.pairs}
\title{ A modified \code{"pairs"} plot with all axes haveing the same range. }
\description{
 Creates all pairwise scatter plots.
}
\usage{
merror.pairs(df,labels=names(df))
}
\arguments{
  \item{df}{ n (no. of items) x N (no. of methods) matrix or data.frame containing the measurements. N must be >= 3 and n > N. }
  \item{labels}{Provide labels for each device down the diagnoal of the pairs plot.}
}
\details{Creates all pairwise scatter plots with the same range for all axes and adds
 the diagonal line denote the "line of equality" or "no bias".).
}
\value{
Produces a scatter plot with the calibration curve and titles that include the
calibration equation and the scale-bias adjusted imprecision standard deviations.
}
\author{ Richard A. Bilonick }

\seealso{ \code{\link{panel.merror} }}
\examples{

data(pm2.5)

# All pairwise plots after square root transformation to Normality
merror.pairs(sqrt(pm2.5))

}

\keyword{ htest }


