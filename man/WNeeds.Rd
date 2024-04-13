\name{WNeeds}
\alias{WNeeds}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Crop water requirement
}
\description{
Calculates the water requirement [m3/ha] of herbaceous or arboreal crops}
\usage{
WNeeds(x, frvol, R)
}

\arguments{
  \item{x}{ Bivariate series of monthly cumulative rainfall [mm] and average monthly temperatures [degree]
}
  \item{frvol}{
Volume fraction of the soil. It is 0.10 for sandy soil, 0.20 fpr loamy soil, 0.18 for clayey soil, 0.13 for medium-textured soil
}
  \item{R}{Length of plant roots [cm] --- see FAO-24 Mannini reworked, maximum extraction depth
}
}

\references{

United States Department of Agricolture (USDA - SCS). IRRIGATION - National Engineering handbook.

Moisello U.  "Idrologia Tecnica" La Goliardica Pavese.

Genovesi R., Bottau D. "L'importanza della falda nell' alimentazione idrica delle colture nella pianura emiliano-romagnola."

Regione Campania - Assessorato Agricoltura - Settore S.I.R.C.A. La tessitura del suolo (foglio divulgativo novembre - dicembre 2002)

Fao irrigation and dreinage paper N. 24 - Crop water requirement, Food and agriculture organization of the united nations ROME, rivisited 1977

Grimaldi, S. Tallerini, C., Serinaldi, F.,  "Modelli multivariati lineari per la generazione di serie di precipitazioni giornaliere", Giornata di Studio: Metodi Statistici e Matematici per l'Analisi delle Serie Idrologiche, Napoli, maggio 2004

}
\author{ Corrado Tallerini
}

\examples{
## data(Pistoia)
## A1=WNeeds(Pistoia,0.13,60)
## edit(A1)
}
