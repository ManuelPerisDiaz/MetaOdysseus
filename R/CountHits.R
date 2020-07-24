#' Count Hit ratio
#'
#' @param TheoreticalMass Theoretical massess
#' @param maldiTOF Output obtained from \code{\link{annotationMALDI_Bot}}.
#' @examples
#' @return Count Hit ratio for PMF MALDI-MS.
#' @export
CountHits<-function(TheoreticalMass,maldiTOF){

total<-as.numeric()
for (i in 1:length(TheoreticalMass)){
 length<-length(TheoreticalMass[[i]][,1])
 total<-sum(as.numeric(length),total)
 }

NumberIonsmatch<-length(maldiTOF$`m/z`)
NumberIonsSubmitted<-total

HitsRatio<-NumberIonsmatch/NumberIonsSubmitted
return(HitsRatio)
}

