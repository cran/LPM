\name{lpm}
\alias{lpm}
\title{ Linear Parametric Model }
\description{
  Estimate ARMA and FARMA models, make simulations and ed eventually apply a corrective procedure to rainfall synthetic series. Besides you can remove seasonal components with STL modified method. }
\usage{
lpm(x, p, q, n, smean, svar, outer=0, prob = 0.95, fre = 365, 
fractional = F, Plag = 20, lsign=0.05, n1 = 399, trasfo = F, des = T, rain = F, graph = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Univariate series }
  \item{p}{ AR order }
  \item{q}{ MA order }
  \item{n}{ Number of series to simulate }
  \item{outer}{ Number of outer loops for STL modified method. Default outer = 0  }
  \item{smean, svar}{ Mean and Variance smoothing windows of STL modified method }
  \item{prob}{ Parameter confidence interval. Default prob = 0.95   }
    \item{fre}{ Series frequency. Default fre = 365 (for daily series)  }
  \item{fractional}{ Logical variable: T to apply FARMA model. Default fractional = F }
  \item{Plag}{ Maximum lag of ACF used in the Portmanteau test. Default Plag = 20 }
  \item{lsign}{Portmanteau Test significance level. Default lsign = 0.05  }

  \item{n1}{ Number of parameters of infinite MA model . Default n1 = 399 }
  \item{trasfo}{ Logical variable: T for preventive logarithmical trasformation. Default trasfo = F}
  \item{des}{ Logical variable: T to remove seasonal components. Default des = T }
  \item{rain}{ Logical variable: T to apply the corrective procedure to daily rainfall simulated series. Default rain = F }
  \item{graph}{ Logical variable: T to receive some graphics. Default graph = F }}

\details{ Need integer periodical dataset. Function to complete modelling univariate series. }
\value{ 
  \item{para }{List of estimated parameters}
  \item{res }{Residual series}
  \item{simdes }{List of simulated series without application of corrective procedure }
  \item{sim }{ List of simulated series }
  \item{BIC }{ Bayesian Information criterion index of estimated model}
 
}
\references{ Grimaldi, S., 'Linear parametric models applied on daily hydrological series', Journal of Hydrologic Engineering, Vol.9, No 5, September  2004.

Grimaldi S., F. Napolitano, L. Ubertini, 'A procedure to use linear parametric models for daily rainfall series simulation' 


Brockwell, P.J and Davis, R.A. (1990) Time Series: Theory and Methods 2nd edition, Springer, NY.

	Hipel, K.W. and McLeod, A.I., (1994) Time Series Modelling of Water Resources and Enviromental Systems, Reading, UK.

	


}
\author{ Salvatore Grimaldi }
\note{ Portmonteau test and BIC index are displaied during application. Portmonteau Test is positive if Q < chi square }

 

\seealso{  \code{\link{rain.adapt}} }
\examples{
##---  lpm(series.runoff,1,1,0,30,30,fractional=T,trasfo=T)
##--   Apply a FARMA(1,d,1) model to series.runoff after e preventive 
##     logarithmical trasformation and deseasonalization with smoothing 30. 

}

