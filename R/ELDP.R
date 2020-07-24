#' Calculates ELDP
#'
#' @param sequence Protein sequence.
#' @param maldiTOF Output obtained from \code{\link{annotationMALDI_Bot}}.
#' @examples
#' @return Calculates ELDP for PMF MALDI-MS.
#' @export


ELDP<-function(maldiTOF,sequence){


  ff<-c(maldiTOF$Seq)
  mt3clv<-OrgMassSpecR::Digest(sequence, missed=2,enzym = "trypsin",IAA=FALSE)


  matched2<-c()
  Results=NULL
  for (i in 1:length(ff)){
    matched<-c()
    a<- match(ff[i],mt3clv[,1] )
    matched <-mt3clv[a,]
    matched2<-rbind(matched2,matched)
  }
  Results<-cbind(maldiTOF,matched2[,2:4])
  Results


  ELDP<-length(which(Results$mc==0))-length(which(Results$mc>0))

  return(ELDP)
}
