\name{PWN}
\alias{PWN}

\title{ Water requirement
}
\description{ Calculate the monthly irrigation requirement of plants based on  cumulative probability [p] and Daily watering duration of irrigation [h]

}
\usage{
PWN(x1, frvol, R, type = "Mi", p, irr)
}

\arguments{
  \item{x1}{ Bivariate series of monthly cumulative rainfall and average monthly temperatures
}
  \item{frvol}{
Volume fraction of the soil. It is 0.10 for sandy soil, 0.20 fpr loamy soil, 0.18 for clayey soil, 0.13 for medium-textured soil
}
  \item{R}{Length of plant roots [cm] --- see FAO-24 Mannini reworked, maximum extraction depth
}
   \item{type}{ type="Ar" for clayey soil, type="Sb" for sandy soil, type="Mi" for medium-textured soil
}
  \item{p}{ Cumulative probability of plant's water requirement [percent]
}
  \item{irr}{ Daily watering duration of irrigation [h]
}
}

\value{
 \item{Values}{Monthly water requirement values [m3/ha] relating to the cumulative probability indicated (p)}

 \item{Flow}{Irrigation flow [l/s/ha] relating to the daily watering duration (irr) and cumulative probability (p)  }
}
\references{

United States Department of Agricolture (USDA - SCS). IRRIGATION - National Engineering handbook.

Moisello U.  "Idrologia Tecnica" La Goliardica Pavese.

Genovesi R., Bottau D. "L'importanza della falda nell' alimentazione idrica delle colture nella pianura emiliano-romagnola."

Regione Campania - Assessorato Agricoltura - Settore S.I.R.C.A. La tessitura del suolo (foglio divulgativo novembre - dicembre 2002)

Grimaldi, S. Tallerini, C., Serinaldi, F.,  "Modelli multivariati lineari per la generazione di serie di precipitazioni giornaliere", Giornata di Studio: Metodi Statistici e Matematici per l'Analisi delle Serie Idrologiche, Napoli, maggio 2004
}
\author{Corrado Tallerini
}

\examples{
##---- data(Pistoia)
##---- PWN(Pistoia,0.13,40,type ="Mi",75,16)
##---- Calculate the monthly irrigation requirement of a plant (Length of plant roots 40 cm in
##---- a medium-textured soil) based on a 75% cumulative probability and 16 hours daily irrigation
}
