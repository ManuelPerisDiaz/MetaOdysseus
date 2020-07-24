#' Scoring function for peptide mass fingerprint MALDI-MS experiment.
#'
#' @param sequence Protein sequence.
#' @param maldiTOF Output obtained from \code{\link{annotationMALDI_Bot}}.
#' @param TheoreticalMass Output obtained from \code{\link{Protein.to.Peptide}}.
#' @param Total Number of amino acids in the protein sequence.
#' @param Mw Molecular weight of the protein assayed.
#' @param write Write output with PMF score.
#' @param enzym Enzyme used for proteolysis.
#' @param missed Number of missed proteolytic cleavages
#' @param save.name Name to save the output results.
#' @examples
#' @return Annotation score for PMF MALDI-MS.
#' @export


ScoresPMF<-function(sequence,maldiTOF,TheoreticalMass,Total,Mw,write=TRUE,save.name,missed,enzym = "trypsin"){

  Sequence.Coverage<-SequenceCounting(sequence,maldiTOF,missed=missed)
  if (length(Sequence.Coverage)==0){Sequence.Coverage<-0
 }
  HitRatio<-CountHits(TheoreticalMass,maldiTOF)
  MassCoverage<-(Sequence.Coverage/100)*Mw
  ELDP <- ELDP(maldiTOF,sequence)
  PMFscore<-(HitRatio*100)+MassCoverage
  if ((Sequence.Coverage)==0){
  PMFscore<- -2 }
  msbcc<-cbind(Sequence.Coverage, HitRatio,MassCoverage, ELDP,PMFscore)

  colnames(msbcc) <- c("Seq. Coverage","Hit Ratio","Mass Coverage","ELDP","PMFscore")

  if (write==TRUE){

    write.csv2(msbcc,paste0("ScoresPMF",save.name,".csv"))
  }

  return(msbcc)

}
