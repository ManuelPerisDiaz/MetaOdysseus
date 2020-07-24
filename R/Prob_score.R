#' Calculate score probability threshold
#'
#' @param maldiTOF Output obtained from \code{\link{annotationMALDI_Bot}}.
#' @param TheoreticalMass Output obtained from \code{\link{Protein.to.Peptide}}.
#' @param mzw Mass error window.
#' @param t2 Left mass
#' @param t1 Right mass
#' @examples
#' @return Calculate score probability threshold for PMF MALDI-MS.
#' @export

Prob_score<-function(maldiTOF,TheoreticalMass,mzw=2,t2=3500,t1=500){
total<-as.numeric()
for (i in 1:length(TheoreticalMass)){
  length<-length(TheoreticalMass[[i]][,1])
  total<-sum(as.numeric(length),total)
}

NumberIonsmatch<-length(maldiTOF$`m/z`)
NumberIonsSubmitted<-total

probabi<-binom.test(x=NumberIonsmatch, n=NumberIonsSubmitted, p=(2*mzw/(t2-t1)))

Prob<-as.numeric(probabi$estimate)
score_prob<- (-log(Prob))


score_prob_threshold<- score_prob
return(score_prob_threshold)

}
