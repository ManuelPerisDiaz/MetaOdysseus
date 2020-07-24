#' Validation top-down MS/MS results
#'
#' @param ev50_true results from topdown.R function
#' @param sequence Protein sequence
#' @param mzw Mass error window.
#' @param Deconv_CysMpro Deconvolved ESI-MS
#' @param Npermu Number of permutations
#' @param intervalmz interval m/z
#' @param Nmetal Number of metal ions
#' @param Metal Metal
#' @examples
#' @return Validation top-down MS/MS results
#' @export
#'
Validation_topdown<-function(ev50_true,sequence,mzw=2,Deconv_CysMpro,

                             intervalmz=3000,Npermu=1000,
                             Nmetal=7,Metal=c("Zn")


){



  theoryMultiple<-MS.MS.Multiple(sequence,native=TRUE,Metal=Metal,Nmetal=Nmetal,Adduct=NULL)

  threshold_b<-sum(ev50_true$score)
  threshold<-  sum(ev50_true$expt_intensity)

  total<-as.numeric()
  for (i in 1:length(  theoryMultiple[[1]])){
    length<-length(  theoryMultiple[[1]][[i]][,1])
    total<-sum(as.numeric(length),total)

  }

  NumberIonsmatch<-length(ev50_true$expt_mz)
  NumberIonsSubmitted<-total
  probabi<-binom.test(x=NumberIonsmatch, n=NumberIonsSubmitted, p=(2*mzw)/intervalmz)
  Prob<-as.numeric(probabi$estimate)
  score_prob_threshold<- (-log(Prob))


  peptideVector <- strsplit(sequence, split = "")[[1]]
  Npermu=Npermu

  Permuted_fasta_lista<-c()
  for(i in 1:Npermu){
    Permuted_fasta <- permute(peptideVector)
    Permuted_fasta_lista[i]<-list(Permuted_fasta)

  }
  pmfscore_totalb<-c()
  pmfscore_total<-c()
  count=0
  count_b=0

  for(i in 1:length(Permuted_fasta_lista)){
    pmfscore<-c()
    pmfscoreb<-c()
    sequence<-paste(Permuted_fasta_lista[[i]],collapse="")

    theoryMultiple<-MS.MS.Multiple(sequence,native=TRUE,Metal=Metal,Nmetal=Nmetal,Adduct=NULL)

    final<-target.MSMS.Multiple(permutation=FALSE,theoryMultiple,Deconv_CysMpro,t=mzw,b=0,native=TRUE,plot= FALSE)


    results<-final
    resb<-results[!duplicated(results$expt_mz),]

    total<-as.numeric()
    for (i in 1:length(  theoryMultiple[[1]])){
      length<-length(  theoryMultiple[[1]][[i]][,1])
      total<-sum(as.numeric(length),total)

    }


    NumberIonsmatch<-length(resb$expt_mz)
    NumberIonsSubmitted<-total

    probabi<-binom.test(x=NumberIonsmatch, n=NumberIonsSubmitted,  p=(2*mzw)/intervalmz)

    Prob<-as.numeric(probabi$estimate)
    score_prob<- (-log(Prob))

    if (is.null(nrow(results))){
      pmfscore<- -20
      pmfscoreb<- 0
    }else{

      pmfscoreb<-score_prob

      pmfscore<-sum(resb$score)


    }

    if ( abs(pmfscore)>  threshold & abs(pmfscoreb)> threshold_b){

      count=count+1

    }


    pmfscore_total<-c(pmfscore_total,(pmfscore))


    if ( score_prob>  score_prob_threshold){

      count_b=count_b+1

    }


    pmfscore_totalb<-c(pmfscore_totalb,(pmfscoreb))


  }

  return(pmfscore_totalb)

}
