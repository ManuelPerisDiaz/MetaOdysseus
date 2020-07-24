#' Determine confidence interval by bootstraping
#'
#' @param data Data
#' @param TheoreticalMass Output obtained from \code{\link{Protein.to.Peptide}}.
#' @param R NUmber of bootstraps
#' @param Mw Molecular weight
#' @param ScorePMF Score value
#' @param sequence protein sequence
#' @param missed Number of missed proteolytic cleavages
#' @examples
#' @return Determine confidence interval by bootstraping for PMF MALDI-MS.
#' @export

PMF_bootstrapping<-function(data,R=999,sequence,TheoreticalMass,Mw,ScorePMF,missed){


rsq<-function(data,indices){
  d<-data[indices,]
  Sequence.Coverage<-SequenceCounting(sequence,d,missed)
  HitRatio<-CountHits(TheoreticalMass,d)
  MassCoverage<-(Sequence.Coverage/100)*Mw
  PMFscore<-(HitRatio*100)+MassCoverage

  pmfscore<-as.numeric(ScorePMF[1,5])
  return(PMFscore)

}

bootcorr <-boot(data=data, statistic = rsq, R=999)

return(bootcorr)

}
