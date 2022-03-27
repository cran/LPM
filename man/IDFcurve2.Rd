\name{IDFcurve2}
\alias{IDFcurve2}
\title{ Intensity duration frequency curve for maximum annual rainfall series of different duration}
\description{
Estimate IDF curve fitting a [mm/h], m ,n, b[h] parameters of maximum annual rainfall series }
\usage{
IDFcurve2(rain, tc, stvalue1 = 1, stvalue2 = 1, t, Tr = 200, MP = F, Trplot = F)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rain}{ Observed Maximum annual rainfall series [mm] of increasing duration  }
   \item{tc}{ Time of concentration of Basin [h] , maybe h(t) and i(t) duration must be calculated }
  \item{stvalue1,stvalue2}{ Starting values of estimation algorithm. Deault stvalue1=1, stvalue2=1 }
  \item{t}{observed rainfall series duration [h] example t=c(1,3,6,12,24) for durations 1,3,6,12,24 hours}
  \item{Tr}{ Return period [y]. Default Tr=200}
   \item{MP}{logical: TRUE for 3 parameters formula i=a/(b+t)^m , FALSE for 2 parameters formula  i=a*t^(n-1), Default MP=False }
   \item{Trplot}{logical: TRUE for plotting Tr values of a(Tr) parameter. Default Trplot=False }
}
\details{
  Estimate parameters of Intensity Duration Frequency curves for maximum annual rainfall series of different duration }
\value{
 \item{par}{List of estimated parameters: a(Tr), m, b, h(t) [mm], i(t) [mm/h], Offset of least squares optimizer}
 \item{I}{I(t) curve scattered point matrix [mm/h]  }
  \item{Curve }{IDF curve scattered point matrix [mm/h] }
  }
\author{ Corrado Tallerini }
\note{a(Tr) is defined by Gumbel distribution.  }
\seealso{  \code{\link{IDFcurve}}}

\examples{
##     ris=IDFcurve2(milano, 1, stvalue1 = 1, stvalue2 = 1,
##     t=c(0.25,0.5,0.75,1,1.25,1.5,2,2.5,3,4,6), Tr = 200, MP=F)
## --> 2 parameters IDF curve estimation of annual maximum rainfall
##     series recorded in Palazzo Marino - Milan (Italy)
##     ris=IDFcurve2(milano, 1, stvalue1 = 1, stvalue2 = 1,
##     t=c(0.25,0.5,0.75,1,1.25,1.5,2,2.5,3,4,6), Tr = 200, MP=T)
## --> 3 parameters IDF curve estimation of annual maximum rainfall
##     series recorded in Palazzo Marino - Milan (Italy)
## --> It's obvious the best performance of the 3 parameters formula}
