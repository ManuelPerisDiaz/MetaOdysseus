#' Protein MS/MS library construction.
#'
#' Calculates y- and b-serie ions for the permuted protein sequences.
#'
#' @param listSeq Output obtain from \code{\link{ProcessMALDIMS}}.
#' @param native If TRUE the library contains metal-protein complexes
#' @param Metal Metal ion considered (e.g. Zn, Cd)
#' @param Nmetal Number of metal ions
#' @param Adduct Adduct
#' @examples
#' @return MS/MS ions corresponding to y- and b-serie from permuted Cys-protein sequences.
#' @export
#'
MS.MS.Multiple<-function(listSeq,native=FALSE,Metal,Nmetal,Adduct){
results<-c()
for (i in 1:length(listSeq)){
theory<-c()
if(native==TRUE){

  Nmetalseq<-seq(0,Nmetal)
  for(z in Nmetalseq){

  theory<-MSMS2(listSeq[[i]],Metal=Metal,Nmetal=Nmetalseq[z+1],Adduct=Adduct)
  results[[i]][z+1]<-list(theory)
  }

}else{
theory<-MSMS(listSeq[[i]])
results[i]<-list(theory)
}}
return(results)
}
