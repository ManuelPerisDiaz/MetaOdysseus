#' Pre-process MALDI-MS spectra.
#'
#' Pre-processing and peak picking for ESI-MS spectra.
#'
#' @param file Experimental MALDI-MS file.
#' @param t1 Select lower m/z limit.
#' @param t2 Select higher m/z limit.
#' @param norm Normalize MS spectra. Default to TRUE.
#' @param comb Select pre-processing method.
#' @param SNR.Th Signal-to-noise ratio. Default to 3.
#' @param amp.Th Parameter controlling the peak detection.
#' @param plotDetection Plot detected peaks after peak picking. Default to FALSE.
#' @param plotProcess Plot processed MS spectra. Default to FALSE.
#' @param plotOriginal Plot raw MS spectra. Default to TRUE.
#' @param dir Directory where to save plots. Default to C:.
#' @param xy MS spectra in xy file?. Default to FALSE.

#' @examples
#' @return Processed MALDI-MS spectra.
#' @export

ProcessMALDIMS<-function(file,
t1=3000,t2=8000 ,comb="comb3",SNR.Th = SNR.Th,amp.Th=amp.Th, xy=FALSE,plotProcess=FALSE, norm=TRUE,plotDetection=FALSE,plotOriginal=TRUE,dir=getwd()){

  if (xy==FALSE){
  filex<-openMSfile(file)
  spec1<-spectra(filex,1)
  peak<-spec1



  ESIprocessed<-ProcessESI(file=file,amp.Th=amp.Th,t1=t1,t2=t2,xy=xy,norm=norm,plotProcess=plotProcess,plotOriginal=plotOriginal,plotDetection = plotDetection,comb=comb,SNR.Th=SNR.Th,dir=dir)

  ESIprocessed<-data.frame(ESIprocessed)
  MatrixMASS.INT<-data.frame(ESIprocessed$mz ,ESIprocessed$V2)
  original<-data.frame(peak[,1],peak[,2])
  return(list("IN.Mz"=MatrixMASS.INT,"OR"=original))

  }else{
    ESIprocessed<-ProcessESI(file=file,amp.Th=amp.Th,t1=t1,t2=t2,xy=xy,norm=norm,plotProcess=plotProcess,plotOriginal=plotOriginal,plotDetection = plotDetection,comb=comb,SNR.Th=SNR.Th,dir=dir)

    ESIprocessed<-data.frame(ESIprocessed)
    MatrixMASS.INT<-data.frame(ESIprocessed$mz ,ESIprocessed$V2)
    return(list("IN.Mz"=MatrixMASS.INT))







  }}


