#' Annotate MALDI-MS spectra
#'
#' Annotate MALDI-MS spectra
#'
#' @param ProcessedMALDIMS Output obtain from \code{\link{ProcessMALDIMS}}.
#' @param database Output obtained from \code{\link{Protein.to.List}}.
#' @param mzw Mass-to-charge window (Da) for the annotation. Default to 1.
#' @param write Write an output table with results. Default to  FALSE.
#' @param orig Plot raw MS spectra. Default to TRUE.
#' @param unique Control annotation. Default to TRUE.
#' @param filterIntensity Filtering by relative intensity. Default to TRUE.
#' @param In Value for filtering by  relative intensity. Default to 10.
#' @param dir Directory where to save plots. Default to C:.
#' @examples
#' \dontrun{
#' MT2a<-("MDPNCSCAAGDSCTCAGSCKCKECKCTSCKKSCCSCCPVGCAKCAQGCICKGASDKCSCCA")
#'TheoreticalMass<-Protein.to.List(MT2a, ET= TRUE)
#'ProcessMALDIMS<-ProcessMALDIMS(file,t1=6000,t2=7000,plotProcess=TRUE,plotDetection=TRUE,plotOriginal = TRUE)
#' maldiTOF<-Target(ProcessMALDIMS,TheoreticalMass, mzw = 2,orig=TRUE,write=FALSE)
#'
#' }
#' @return Peaks annotated in the MALDI-MS spectra.
#' @export

Target<- function (ProcessedMALDIMS,database, mzw = 1, write=FALSE,orig=TRUE,unique=TRUE,filterIntensity=TRUE,In=10,dir="C:" ){
  exp=NULL
  or=NULL
  y=NULL
  yx=NULL

exp<-as.data.frame(ProcessedMALDIMS$IN.Mz)
or<-as.data.frame(ProcessedMALDIMS$OR)
y<-exp[,2]
yx = (y-min(y))/(max(y)-min(y))*100


 exp<-cbind(exp,ppmdif="ppmdif")
 exp<-cbind(exp,N.Modifications="N.Modifications")
 exp<-cbind(exp,yx)
  f<-mzw
  ions1bb<-data.frame()
  msbc<-data.frame()
  b<-data.frame()
  msb<-data.frame()
  ms <- data.frame()
  for (i in 1:nrow(database)){

    from<-2
    to<-as.numeric(length(database[i,]))
    fullmsset <- database[i, from:to]
    options(digits=6)
    fullms2 <-  as.numeric(fullmsset[!is.na(fullmsset)])
    name<- database[i, 1]

    msb <- exp[which(exp[, 1] >  fullms2 - f & exp[, 1] <  fullms2 + f ),]

      if (nrow(msb) > 0){
        ppmdif<-((as.numeric(msb[1])- fullms2))
        msb[,3]<- ppmdif

        msb[,4]<-name
    }
          msbc<-rbind(msb,msbc)
    msbc<-as.data.frame(msbc)
  }

    if(nrow(msbc) > 0){
     if(write==TRUE){
       colnames(msbc) <- c("m/z","Maximum Intensity","Error(Da)","N.Modifications", "Rel.Intensity")

      write.csv2(msbc,paste0("MALDIannotated",".csv"))}
  }else{return(print("No peaks matched"))

  }


  colnames(msbc) <- c("m/z","Maximum Intensity","Error(Da)","N.Modifications", "Rel.Intensity")
if(unique==TRUE){



  aa<-msbc[order(msbc$`Maximum Intensity`,-abs(msbc$`Maximum Intensity`)),]

  aaa<-aa[!duplicated(aa$N.Modifications),]

  msbc<-aaa[order(aaa$N.Modifications),]

}




   subsettedInten<-subset(msbc, msbc[,5]>In)
   FiltretedQC1<-cbind(  subsettedInten[,1],subsettedInten[,5])
   plot_detection(FiltretedQC1,orig=FALSE,do.plot = FALSE,QC1=TRUE,QC2=FALSE,norm=TRUE,dir=dir)
   FiltretedQC2<-cbind(subsettedInten[,3],subsettedInten[,5])
   plot_detection(FiltretedQC2,orig=FALSE,do.plot = FALSE,QC1=FALSE,QC2=TRUE,norm=TRUE,dir=dir)


    return(subsettedInten)




  }





