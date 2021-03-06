% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Target.R
\name{Target}
\alias{Target}
\title{Annotate MALDI-MS spectra}
\usage{
Target(
  ProcessedMALDIMS,
  database,
  mzw = 1,
  write = FALSE,
  orig = TRUE,
  unique = TRUE,
  filterIntensity = TRUE,
  In = 10,
  dir = "C:"
)
}
\arguments{
\item{ProcessedMALDIMS}{Output obtain from \code{\link{ProcessMALDIMS}}.}

\item{database}{Output obtained from \code{\link{Protein.to.List}}.}

\item{mzw}{Mass-to-charge window (Da) for the annotation. Default to 1.}

\item{write}{Write an output table with results. Default to  FALSE.}

\item{orig}{Plot raw MS spectra. Default to TRUE.}

\item{unique}{Control annotation. Default to TRUE.}

\item{filterIntensity}{Filtering by relative intensity. Default to TRUE.}

\item{In}{Value for filtering by  relative intensity. Default to 10.}

\item{dir}{Directory where to save plots. Default to C:.}
}
\value{
Peaks annotated in the MALDI-MS spectra.
}
\description{
Annotate MALDI-MS spectra
}
\examples{
\dontrun{
MT2a<-("MDPNCSCAAGDSCTCAGSCKCKECKCTSCKKSCCSCCPVGCAKCAQGCICKGASDKCSCCA")
TheoreticalMass<-Protein.to.List(MT2a, ET= TRUE)
ProcessMALDIMS<-ProcessMALDIMS(file,t1=6000,t2=7000,plotProcess=TRUE,plotDetection=TRUE,plotOriginal = TRUE)
maldiTOF<-Target(ProcessMALDIMS,TheoreticalMass, mzw = 2,orig=TRUE,write=FALSE)

}
}
