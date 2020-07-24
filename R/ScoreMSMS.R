#'  Scoring function for MS/MS spectra.
#'
#' @param theoryMultiple Theoretical MS/MS library.
#' @param final Output Theoretical MS/MS library.
#' @param exp Experimental MS/MS file.
#' @param tol Mass error
#' @param exp Experimental MS/MS file.
#' @param Score Select the type of score using for finding best prediction following the next order: These are number of ions matched, percentage peaks matched and percentage of ions matched.
#' @param write.csv Write output.
#' @param write.csvb Write output.
#' @examples
#' @return Annotation score for MS/MS spectra
#' @export

ScoreMSMS<-function(theoryMultiple,final,tol,exp,Score=4,write.csv=TRUE,write.csvb=TRUE){


  Intensityionsmached_global<-c()
  Percentofionsmached_global<-c()
  Percentofpeaksmached_global<-c()
  NumberIonsmatch_global<-c()
  MeanError_global<-c()
  Score_global<-c()
  Results<-matrix(,ncol=7,nrow=length(theoryMultiple))

for (i in 1:length(theoryMultiple)){
the<-theoryMultiple[[i]]
target<-final[[i]]
if(!is.null(target)){

MeanError<-mean(target$error)
NumberIonsmatch<-length(target$ms2mz)
Percentofpeaksmached<-(length(target$ms2mz)/length(the$ms2mz))*100
Percentofionsmached<-(length(target$ms2mz)/length(exp[,1]))*100
Intensityionsmached<-((sum(target$expt_intensity))/(sum(exp[,2])))*100


NumberIonsSubmitted<-length(the$ms2mz)
mzw = tol
t1=min(exp[,1])
t2=max(exp[,1])
probabi<-binom.test(x=NumberIonsmatch, n=NumberIonsSubmitted, p=2*mzw/(t2-t1),alternative="two.sided")
Prob<-as.numeric(probabi$p.value)
score_prob<- (-log(Prob))


Results[i,1]<-((target$ms1seq))[2]
Results[i,2]<-as.numeric(MeanError)
Results[i,3]<-as.numeric(NumberIonsmatch)
Results[i,4]<-as.numeric(Percentofpeaksmached)
Results[i,5]<-as.numeric(Percentofionsmached)
Results[i,6]<-as.numeric(Intensityionsmached)
Results[i,7]<-score_prob

MeanError_global<-c(MeanError_global,MeanError)
Intensityionsmached_global<-c(Intensityionsmached_global,Intensityionsmached)
Percentofionsmached_global<-c(Percentofionsmached_global,Percentofionsmached)
Percentofpeaksmached_global<-c(Percentofpeaksmached_global,Percentofpeaksmached)
NumberIonsmatch_global<-c(NumberIonsmatch_global,NumberIonsmatch)

Score_global<-c(Score_global,score_prob)

}
}
  colnames(Results) <- c("Sequence","Mean Error (Da)","N. ions matched","Per. peaks matched","Per. ions matched","Int. ions matched","Score")
  Results<-data.frame(Results)

  Results<- Results[complete.cases(Results), ]

  if (write.csv==TRUE){
    write.csv2(Results,file="Score.csv")
  }

  if (Score==1){
  i<-which.max(Results[,4])
  SeqRe<-Results[i,]
  }
  if (Score==2){
    i<-which.max(Results[,5])
    SeqRe<-Results[i,]
  }
  if (Score==3){
    i<-which.max(Results[,3])
    SeqRe<-Results[i,]
  }

  if (Score==4){
    i<-which.max(Results[,7])
    SeqRe<-Results[i,]
  }
  if (write.csvb==TRUE){
    write.csv2(SeqRe,file="SeqRe.csv")
  }



  return(SeqRe)

  }
