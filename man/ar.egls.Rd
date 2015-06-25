\name{ar.egls}
\alias{ar.egls}
\title{ Subset Autoregressive Model }
\description{
Estimate VAR(p) model fixing some parameter values to zero}
\usage{
ar.egls(x, R, order.max , na.action = na.fail, series = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Univariate or multivariate series with nil mean }
  \item{R}{ Matrices of parameters selection }
  \item{order.max}{ Model order }
  \item{na.action}{ Function to be called to handle missing values }
  \item{series}{ Names for the series.  Defaults to 'deparse(substitute(x))'}
  \item{\dots}{ See ar.ols }
}
\details{
  R matrix is a list of p matrices, with p the autoregressive order. In R value '1' allows parameter estimation, '0' fix the parameter value to zero. }
\value{ See ar.ols
  
}
\references{ Grimaldi  S. , Serinaldi F. & Tallerini C. (2004)  'Multivariate linear parametric models applied to daily rainfall time series'  Mediterranean Storms, 6rd EGU Plinius Conference held in Mediterranean Sea, Italy, October 2004
 
Lutkepohl, H. (1993) Introduction to Multiple Time Series
     Analysis 2nd Edition ._ Springer Verlag, NY}
\author{ Corrado Tallerini }
\note{ Function is created modifing ar.ols by Adrian Trapletti and Brian Ripley }


\examples{
##	S1=matrix(0,3,3)
##	S1[1,1]=1
##	S1[1,2]=1
##	S=list()
##	S[[1]]=S1
##	S[[2]]=S1
##	ar.egls(series.rainfall[,1:3],S,order.max=2)
## --> Apply a Subset VAR(2) model restricted to 4 parameters (position (1,1)
## --> and (1,2) in both matrices) to first 3 series of series.rainfall 
## --> dataset  }

