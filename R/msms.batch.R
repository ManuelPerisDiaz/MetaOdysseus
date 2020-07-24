#' Produce a decoy database
#'
#' @param sequence protein sequence
#' @param missed Control enzymatic digestion
#' @param tol mass error tolerance
#' @param mzML If TRUE, the file input is in mzML file format
#' @param label Name to save
#' @param Score Score
#' @param plot If TRUE, plot results
#' @param write.csv IF TRUE, expor results in csv
#' @param write.csvb IF TRUE, expor results in csv
#' @examples
#' @return Decoy database
#' @export
#'
msms.batch<-function(sequence,missed,tol,mzML=TRUE,label="Decoy",Score,
                     plot=FALSE,write.csv=FALSE,write.csvb=TRUE){



  if (mzML==TRUE){
    folder<-dir()[grep(".mzML",dir())]


  }else{

    folder<-dir()[grep(".txt",dir())]

  }


  matches_f <- vector("list")
  parent_f<-c()
  TheoreticalMass<-Protein.to.Peptide(sequence,missed=missed)


  for (z in 1:length(folder)){
    exp<-c()
    parent<-c()
    Res<-c()

    if (mzML==TRUE){

      expb<-openMSfile(folder[z])
      exp<-spectra(expb,1)
      parent<- header(expb,1)$precursorMZ

    }else{

      exp<-read.table(folder[z])
      parent<-as.numeric(regmatches(folder[z], gregexpr("[[:digit:]]+", folder[z]))[[1]][1]
      )
    }



    dir.create(paste(parent[1],label,sep=""),showWarnings = FALSE)
    options(warn=-1)


    ScoreMSMSre<-ms.ms_bottom(sequence,exp=exp,missed=missed,parent=parent,tol=tol,
                       Score=Score,plot=plot,write.csv=write.csv,write.csvb=write.csvb)

    if (!is.na(ScoreMSMSre)){
    ScoreMSMSre<-cbind(ScoreMSMSre,parent)

    write.csv2(ScoreMSMSre,paste(parent[1],label,".csv",sep=""))


    new.folder <-paste(parent[1],label,sep="")
    file.copy(folder[z], new.folder)
    file.copy("Score.csv", new.folder)
    file.copy(paste(parent[1],label,".csv",sep=""), new.folder)

    parent_f<-c(parent_f, parent)

    matches_f[[z]] <-data.frame(sequence=as.character(ScoreMSMSre$Sequence),Mean.Error..Da.=as.numeric(as.character(ScoreMSMSre$Mean.Error..Da.)), N..ions.matched=as.numeric(as.character(ScoreMSMSre$N..ions.matched)),

                                Per..peaks.matched=as.numeric(as.character(ScoreMSMSre$Per..peaks.matched)), Per..ions.matched=as.numeric(as.character(ScoreMSMSre$Per..ions.matched)),
                                Int..ions.matched=as.numeric(as.character(ScoreMSMSre$Int..ions.matched)),Score= as.numeric(as.character(ScoreMSMSre$Score)),parent= as.numeric(as.character(ScoreMSMSre$parent)))


    matches_f[[z]] <-data.frame(sequence=ScoreMSMSre$Sequence,Mean.Error..Da.=ScoreMSMSre$Mean.Error..Da., N..ions.matched=ScoreMSMSre$N..ions.matched,

                                Per..peaks.matched=ScoreMSMSre$Per..peaks.matched, Per..ions.matched=ScoreMSMSre$Per..ions.matched,
                                Int..ions.matched=ScoreMSMSre$Int..ions.matched,Score= ScoreMSMSre$Score,parent= as.factor(ScoreMSMSre$parent))




    if(plot==TRUE){

      image<-paste0("theoryMultiple",1,".tiff",sep="")
      file.copy(image, new.folder)
    }

  }}
  return( matches_f)
}
