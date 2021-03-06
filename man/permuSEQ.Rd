% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/permuSEQ.R
\name{permuSEQ}
\alias{permuSEQ}
\title{Permuted protein sequence.}
\usage{
permuSEQ(sequence, NModif = 3)
}
\arguments{
\item{sequence}{Protein sequence.}

\item{NModif}{Number of total modifications expected in the Cys-protein residues.}
}
\value{
Permuted Cys-protein sequences.
}
\description{
Calculates permuted protein sequences according to number of Cys residues present.
}
\examples{
\dontrun{
sequence<-c("KSCCSCCPAECEK")
listSeq<-permuSEQ(sequence,NModif=4)
}
}
