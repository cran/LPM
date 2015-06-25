\name{mlpm}
\alias{mlpm}

\title{ Multivariate Linear Parametric Model }
\description{
 Multivariate modelling using VAR(p) and SVAR(p) different estimation methods, simulation, daily rainfall simulated series correction and deseasonalization are performed

  }
\usage{
mlpm(x, p, prob, nsim, smean, svar, fre = 365, outer = 0,plot = F, 
rain = T, over = T, estimate = "ols", CCFlag = 20, Plag = 20, lsign = 0.05, des = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Multivariate series }
  \item{p}{ Model order }
  \item{prob}{ Condifidence interval used to fix parameters in SVAR(p) model  }
  \item{nsim}{ Number of series to simulated }
  \item{smean , svar}{ Mean and Variance smoothing windows of STL modified method }
  \item{fre}{ Series frequency. Default fre = 365 }
  \item{outer}{ Outer loops of STL modified method. Default outer = 0 }
\item{plot}{ Logical variable: T to receive some graphics. Default plot = F  }
  \item{rain}{ Logical variable: T to apply rain adaptor to simulated series. Default rain = F  }
  \item{over}{ Logical variable: T to use SVAR(p) model estimated with EGLS method. Need estimate = 'ols'  Default over = T  }
  \item{estimate}{ Define VAR(p) estimation method. 'ols', 'burg', 'yw' (Yule-Walker). Default estimate = 'ols' }
  \item{CCFlag}{ Lag of (Partial) Auto-CrossCorrelation function graphics . Default CCFlag = 20 }
  \item{Plag}{ Maximum lag of A-CCF used in the Portmanteau Test. Default Plag = 20 }
  \item{lsign}{Portmanteau Test significance level. Default lsign = 0.05  }
  \item{des}{ Logical variable: T to remove seasonal components  }
}
\details{
Need integer periodical datasets.
Simulation use Lutkepohl algorithm with a residuals vectorial permutation to obtain innovations.
Parameters selections of EGLS method is defined by t-ratio approach.

}
\value{

  \item{coeff }{List of estimated coefficients matrix}
  \item{coeffstd }{ List of estimated standard deviations coefficients matrix. Only for OLS and EGLS method}

  \item{struct }{ List of 'structure' of SVAR(p) model (1 define position of estimated parameter). Only for EGLS method
}
  \item{res }{Residual series }
  \item{fit }{Output List of ar function}
  \item{aic}{Akaike Information Criterion index }
  \item{Q}{Portmonteau statistic }
  \item{sim}{List of simulated series  }
}


\references{ Grimaldi S., Tallerini C., Serinaldi F. (2004)  'Modelli multivariati lineari per la generazione di serie di precipitazioni giornaliere' Giornata di Studio: Metodi Statistici e Matematici per l'Analisi Idrologiche  Napoli 2004

Grimaldi  S. , Serinaldi F. & Tallerini C. (2004)  'Multivariate linear parametric models applied to daily rainfall time series'  Mediterranean Storms, 6rd EGU Plinius Conference held in Mediterranean Sea, Italy, October 2004

Lutkepohl, H. (1993) Introduction to Multiple Time Series Analysis 2nd edition, Springer-Verlag, Berlin.

Grimaldi, S.,  'Linear parametric models applied on daily hydrological series', Journal of Hydrologic Engineering, Vol. 9, No 5 , September 2004.

Brockwell, P.J and Davis, R.A. (1990) Time Series: Theory and Methods 2nd edition, Springer, NY.

Hipel, K.W. and McLeod, A.I., (1994) Time Series Modelling of Water Resources and Enviromental Systems, Reading, UK.

     Hosking, J.R.M. (1980) 'The Multivariate Portmanteau Statistic' Journal of the American Statistical Association, Vol.75, N.371, 502-608.
}
\author{Corrado Tallerini  }
\note{ Portmonteau test, AIC e SBC index are displaied during application. Portmonteau test is positive if Q < chi square. }


\seealso{  \code{\link{lpm}}, \code{\link{ar.egls}}, \code{\link{rain.adapt}}}
\examples{
##-- Mrain=mlpm(series.rainfall,3,0.95,0,120,120)
##-- Apply a SVAR(3) model with selection probability 95 % to series.rainfall
##-- after preventive deseasonalization with smoothing 120.

}


