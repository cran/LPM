\name{rain.adapt}
\alias{rain.adapt}
\title{ Rainfall Adaptor  }
\description{
  Apply a corrective procedure to daily rainfall series to enforce actual caracteristics.
}
\usage{
rain.adapt(x, a, ser)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Observed series }
  \item{a}{ Univariate series to modify (simulated series) }
  \item{ser}{ Series identification number  }
}
\details{ The no-rain frequency consequentally the total rainfall depth of the observed series are enforced on the synthetic series
}
\value{
Corrected series
}
\references{ Grimaldi S., F. Napolitano, L. Ubertini, 'A procedure to use linear parametric models for daily rainfall series simulation'  }
\author{Salvatore Grimaldi }
\examples{
##   rain=lpm(series.rainfall[,1],1,1,1,120,120)
##   rain.adapt(series.rainfall[,1],rain$sim[[1]],1)
##-- ==>  Apply rain adaptor to a simulated series with a ARMA(1,1) model

}

