#' Validate the scoring function for peptide mass fingerprint MALDI-MS experiment.
#'
#' @param ProcessMALDIMS Processed MALDI-MS spectrum Protein sequence.
#' @param TheoreticalMass Output obtained from \code{\link{Protein.to.Peptide}}.
#' @param nAA Number of amino acids in the protein sequence.
#' @param Mw Molecular weight of the protein assayed.
#' @param mzw Mass error window.
#' @param threshold threshold for scoring
#' @param sequence Protein sequence
#' @param Npermu Number of permutations
#' @param In In for annotationMALDI_Bot function
#' @param missed Number of missed proteolytic cleavages
#' @param prob If TRUE returns the score calculated with binomial function
#' @param score_prob_threshold Threshold for probability
#' @examples
#' @return Annotation score for PMF MALDI-MS.
#' @export


Validate_score<-function(ProcessMALDIMS,TheoreticalMass,mzw,nAA,threshold, Mw,sequence,Npermu=10000,
                         In,prob=FALSE,score_prob_threshold,missed){
peptideVector <- strsplit(sequence, split = "")[[1]]


Permuted_fasta_lista<-c()
for(i in 1:Npermu){
  Permuted_fasta <- permute(peptideVector)
  Permuted_fasta_lista[i]<-list(Permuted_fasta)

}
score_prob_total<-c()
pmfscore_total<-c()
count=0
count_b=0
for(i in 1:length(Permuted_fasta_lista)){
  pmfscore<-c()
  sequence<-paste(Permuted_fasta_lista[[i]],collapse="")

  TheoreticalMass<-Protein.to.Peptide(sequence=sequence)
  maldiTOF<-annotationMALDI_Bot(ProcessedMALDIMS=ProcessMALDIMS,TheoreticalMass=TheoreticalMass,mzw = mzw,plot=FALSE,plotIndi=FALSE,orig=TRUE,write=FALSE,In=In)


  if (is.null(nrow(maldiTOF))){
    pmfscore<- -2

    score_prob<-0
  }else{

    total<-as.numeric()
    for (i in 1:length(TheoreticalMass)){
      length<-length(TheoreticalMass[[i]][,1])
      total<-sum(as.numeric(length),total)
    }

    NumberIonsmatch<-length(maldiTOF$`m/z`)
    NumberIonsSubmitted<-total

    probabi<-binom.test(x=NumberIonsmatch, n=NumberIonsSubmitted, p=2*mzw/3000)

    Prob<-as.numeric(probabi$p.value)
    score_prob<- (-log(Prob))

    CountAAPMFF<-SequenceCounting(seq=sequence,maldiTOF=maldiTOF,missed=missed)

    ScorePMF<-ScoresPMF(sequence=sequence,maldiTOF=maldiTOF,TheoreticalMass=TheoreticalMass,Total=nAA,Mw=Mw,write = FALSE,missed=missed)

    pmfscore<-as.numeric(ScorePMF[1,5])

  }

  if ( pmfscore>  threshold){

    count=count+1

  }

  if ( score_prob>  score_prob_threshold){

    count_b=count_b+1

  }


  pmfscore_total<-c(pmfscore_total,pmfscore)
  score_prob_total<-c(score_prob_total,score_prob)
}

if(prob=="TRUE"){
return(score_prob_total)}else{

  return(  pmfscore_total)
}

}
