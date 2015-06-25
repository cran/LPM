\name{series.rainfall}
\alias{series.rainfall}
\docType{data}
\title{  Daily Rainfall Series }
\description{
Group of 5 daily rainfall series recorded in Tuscany region of Italy during the period 1958-1979. }
\usage{data(series.rainfall)}

\details{
Dataset is created removing lacking years and replacing lacking days with the mean of previous and successive value. Beside 29 february day values are removed   to obtain integer periodical dataset.}
\source{
Rudari, R. 'Predicibilita' del clima europeo ed influenze delle forzanti a scala sinottica su eventi regionali di precipitazione intensa', PDh Thesis 2001}
\examples{
data(series.rainfall)
## maybe str(series.rainfall) ; plot(series.rainfall) ...
}
\keyword{datasets}
