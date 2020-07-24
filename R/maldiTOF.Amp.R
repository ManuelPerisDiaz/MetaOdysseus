#' Extended annotation for Bottom-up MALDI-MS experiments
#'
#' Extended annotation for Bottom-up MALDI-MS experiments including the start, stop and missed cleavage sites from the proteolytic reaction.
#'
#' @param sequence Protein sequence.
#' @param maldiTOF Output obtained from \code{\link{annotationMALDI_Bot}}.
#' @param missed Missed cleavage sites produced by proteolytic enzyme. Default to 2.
#' @param enzyme Proteolytic enzyme used in the experiment. Default to Trypsin.
#' @param write Write results. Default to TRUE.
#' @examples
#' @return Extended annotation for bottom-up MALDI-MS experiment.
#' @export
maldiTOF.Amp<-function(sequence,maldiTOF,missed=2,enzyme="trypsin",write=TRUE){
Results=NULL
ff<-c(maldiTOF$Seq)
mt3clv<-Digest(sequence, missed,enzyme=enzyme ,IAA=FALSE)


matched2<-c()
for (i in 1:length(ff)){
  matched<-c()

  i<- match(ff[i],mt3clv[,1] )
  matched <-mt3clv[i,]
  matched2<-rbind(matched2,matched)
}
Results<-cbind(maldiTOF,matched2[,2:4])

if (write==TRUE){

  write.csv2(Results,"MALDIannotatedAmp2.csv")

}
return(Results)
}
