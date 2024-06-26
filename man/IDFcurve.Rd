\name{IDFcurve}
\alias{IDFcurve}
\title{ Intensity duration frequency curve }
\description{
Estimate IDF curve fitting a [mm/h], m ,n, b[h] parameters }
\usage{
IDFcurve(rain, g, s, tc, stvalue1 = 1, stvalue2 = fre, fre, Tr = 200,
                      MP=F, Trplot=F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rain}{ Observed Univariate rainfall series non cumulative }
  \item{g}{ Maximum bound for cumulative series. For daily series g = 7 is recommended, for hourly series g=24 is racommended }
  \item{s}{ Threshold for defining "event". If  "10", only h > = 10 mm values are considered}
   \item{tc}{ Time of concentration of Basin [h]}
  \item{stvalue1,stvalue2}{ Starting values of estimation algorithm. Deault stvalue1=1, stvalue2=fre }
  \item{fre}{ Series frequency [h]. For daily series fre=24, for hourly series fre=1}
  \item{Tr}{ Return period [y]. Default Tr=200}
  \item{MP}{logical: TRUE for 3 parameters formula i=a/(b+t)^m , FALSE for 2 parameters formula  i=a*t^(n-1), Default MP=False }
   \item{Trplot}{logical: TRUE for plotting Tr values of a(Tr) parameter. Default Trplot=False }
     }
\details{
  Estimate parameters of Intensity Duration Frequency curves }
\value{
 \item{par}{List of estimated parameters: a(tr), m, b, h(t) [mm], i(t) [mm/h], Offset of least squares optimizer}
  \item{Curve }{IDF curve Scattered point matrix [mm/h]  }
  }
\author{ Corrado Tallerini }
\note{a(tr) is defined by Gumbel distribution.  }
\seealso{  \code{\link{IDFcurve2}}}

\examples{
##    data(hourly.rainfall.series)
##    res = IDFcurve(hourly.rainfall.series ,24, 15, 1, fre=1, Tr=200, MP=F)
## --> 2 parameters IDF curve estimation of a hourly rainfall series applying
## --> a Threshold "15 mm"  and Time of concentration t=1 h
##    res = IDFcurve(hourly.rainfall.series ,24, 15, 1, fre=1, Tr=200, MP=T)
## --> 3 parameters IDF curve estimation of a hourly rainfall series applying
## --> a Threshold "15 mm"  and Time of concentration t=1 h
## --> It's obvious the best performance of the 3 parameters formula}
