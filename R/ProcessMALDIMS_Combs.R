ProcessMALDIMS_Combs.R<-function(file,
t1=3000,t2=8000, SNR=2,comb="comb3",plotProcess=TRUE, norm=TRUE,plotDetection=FALSE){


  filex<-openMSfile(file)
  spec1<-spectra(filex,1)
  peak<-spec1

  ESIprocessed<-ProcessESI(file=file,t1=t1,t2=t2,norm=norm,plotProcess=plotProcess,plotDetection = plotDetection,comb=comb,SNR.Th=SNR)

  ESIprocessed<-data.frame(ESIprocessed)
  MatrixMASS.INT<-data.frame(ESIprocessed$mz ,ESIprocessed$V2)
  original<-data.frame(peak[,1],peak[,2])
  return(list("IN.Mz"=MatrixMASS.INT,"OR"=original))

  }


