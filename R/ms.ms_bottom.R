
ms.ms_bottom<-function(sequence,exp,missed=2,parent,tol=2,Score=1,plot,write.csv=TRUE,write.csvb=FALSE){

  TheoreticalMass<-Protein.to.Peptide(sequence,missed=missed)


  theoryMultiple_global_fin<-c()


  for (i in 1:length(TheoreticalMass)){


    nameTry<-names(TheoreticalMass[i])
    database<-TheoreticalMass[[i]]


    theoryMultiple_global<-c( )


    for (i in 1:nrow(database)){
      from<-2
      to<-as.numeric(length(database[i,]))
      fullmsset <- database[i, from:to]
      options(digits=6)
      fullms2 <-  as.numeric(fullmsset[!is.na(fullmsset)])
      NModif<- as.numeric(database[i, 1])

      if (abs(fullms2-parent)<=  tol ){


        if (NModif=="0"){
          listSeq_theo<-nameTry

        }else{
          listSeq_theo<-permuSEQ( nameTry,NModif=NModif)

        }

        if( nchar(listSeq_theo)>1){

          theoryMultiple<-MS.MS.Multiple(listSeq_theo)


          theoryMultiple_global<-c(theoryMultiple_global,theoryMultiple)


        }}


    }


    theoryMultiple_global_fin<-c(theoryMultiple_global_fin,theoryMultiple_global)
  }

  if (is.null(theoryMultiple_global_fin)){
    ScoreMSMS<-NA

  }else{

    final<-target.MSMS.Multiple_B(theoryMultiple=theoryMultiple_global_fin,exp=exp,t=tol,plot=plot,xlim=NULL)

    ScoreMSMS<-ScoreMSMS(theoryMultiple=theoryMultiple_global_fin,final,tol=tol,exp=exp,Score=Score,write.csv=write.csv,write.csvb=write.csvb)


  }
  return(ScoreMSMS)

}


