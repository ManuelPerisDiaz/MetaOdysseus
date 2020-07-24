#' Pre-process ESI-MS spectra.
#'
#' Pre-processing and peak picking for ESI-MS spectra.
#'
#' @param file Experimental MS file.
#' @param amp.Th Parameter controlling the peak detection.
#' @param t1 Select lower m/z limit.
#' @param t2 Select higher m/z limit.
#' @param norm Normalize MS spectra. Default to TRUE.
#' @param QC1 Generate quality control type 1 plot. Default to FALSE.
#' @param QC2 Generate quality control type 2 plot. Default to FALSE.
#' @param comb Select pre-processing method.
#' @param SNR.Th Signal-to-noise ratio. Default to 3.
#' @param plotDetection Plot detected peaks after peak picking. Default to FALSE.
#' @param plotProcess Plot processed MS spectra. Default to FALSE.
#' @param plotOriginal Plot raw MS spectra. Default to TRUE.
#' @param xy MS spectra in xy file?. Default to FALSE.
#' @param dir Directory where to save plots. Default to C:.
#' @return Processed ESI-MS spectra.
#' @examples
#' @export


ProcessESI<-function(file,amp.Th=0.003,t1=300,t2=1800,norm=TRUE,QC1=FALSE, QC2=FALSE,comb="comb3",SNR.Th=2,plotDetection=FALSE,plotProcess=FALSE,plotOriginal=TRUE,xy=FALSE,dir="C:"){

  if (xy==TRUE){
    peak<-data.frame(file)

  }else{
   filex<-openMSfile(file)
  spec1<-spectra(filex,1)
  peak<-spec1
  }

  if (t1 < 0 | t2 < 0 | t1 > t2) {
    stop("No possible values for t1 and t2")
  }


  if (t2 > t1) {
    peakB <- peak[which(peak[, 1] > t1 & peak[, 1] < t2),
                  ]
  }

  smoothDWT  = getFromNamespace("smoothDWT", "MassSpecWavelet")
  sav.gol = getFromNamespace("sav.gol", "MassSpecWavelet")
  if (comb== "comb1") {
    ################ process the spectra, baseline and smooth
    ## comb1 smooth savitsky golay + baseline ALS
    v<-sav.gol(peakB[,2],fl=200)
    ycorr<-baseline.corr(v)
    Resulted<-cbind(peakB[,1],ycorr)
    }

  if (comb== "comb2") {
## comb2  baseline ALS +smooth savitsky golay
ycorr<-baseline.corr(peakB[,2])
v<-sav.gol(ycorr,fl=50)
Resulted<-cbind(peakB[,1],v)
}

  if (comb== "comb3") {
## comb3 smooth DIFSM + baseline ALS
smoothed<-difsm(peakB[,2], lambda = 1e7)
ycorr<-baseline.corr(smoothed)
Resulted<-cbind(peakB[,1],ycorr)
  }

  if (comb== "comb4") {
## comb4  baseline ALS +smooth DIFSM
ycorr<-baseline.corr(peakB[,2])
smoothed<-difsm(ycorr, lambda = 1e7)
Resulted<-cbind(peakB[,1],smoothed)
}
  if (comb== "comb5") {
## comb5  baseline ALS +smooth DWT
ycorr<-baseline.corr(peakB[,2])
smoothed<-smoothDWT(ycorr)
Resulted<-cbind(peakB[,1],smoothed)
  }
if (comb== "none"){

  Resulted<- peakB
}

  dir = dir
  file=paste(dir,"Preprocessed_spectra.csv",sep="")
  write.csv2(Resulted,file)



  if (plotProcess== TRUE) {
    plot_processed(Resulted,norm=norm,dir=dir,csv=FALSE)
    }

  #

################
#### PEAK PICKING

peakInfo <- peakDetectionCWT(Resulted[,2],SNR.Th=SNR.Th,amp.Th = amp.Th)
majorPeakInfo = peakInfo$majorPeakInfo
peakIndex <- majorPeakInfo$peakIndex

mz<-peakB[,1][peakInfo$majorPeakInfo$peakCenterIndex]

DetectedSmoothed<-cbind(mz,peakInfo$majorPeakInfo$peakValue)

if (plotDetection== TRUE) {
  plot_detection(DetectedSmoothed,orig=FALSE,QC2 =FALSE,QC1=FALSE,norm=norm,dir=dir)

}

if (plotOriginal ==TRUE){

  plot_Peak(peakB,norm=norm,dir=dir)
}

return(DetectedSmoothed)

}




