#' Annotate a top-down MS/MS spectra
#'
#' @param experimental Experimental MS/MS spectra
#' @param labelling If TRUE, a CYs labelling reagent is used.
#' @param sequence Protein sequence
#' @param t Mass error window.
#' @param mz1 left mass
#' @param mz2 right mass
#' @param max.charge maximum charge state for deconvolution
#' @param tol mass tolerance
#' @param b Filter by intensity
#' @param SNR.Th SNR threshold
#' @param peak if TRUE peak picking is performed
#' @param plot if TRUE, plot is produced
#' @param Maxcharge maximum charge state for theoretical database
#' @param maxMass maximum mass for deconvolution
#' @param parent parent ion selected for MS/MS
#' @param mzw mass error window to annotate the parent ion
#' @param NbindingMetal Number of metal ions
#' @param Metal Metal
#' @param Adduct If NULL, proton is the adduct considered
#' @examples
#' @return Annotate top-down MS/MS spectra
#' @export
#'
Topdown<-function(experimental,labelling=FALSE, mz1=500,mz2=3000,
                  max.charge,tol,sequence,t, b,SNR.Th,peak=TRUE,plot=FALSE,
                  maxMass,parent,mzw, Metal=c("Zn"), NbindingMetal = 7,
                  Maxcharge=6,Adduct=NULL){

  if(labelling==TRUE){
    ESI.Theoretical_res<-ESI.theoretical.alkylators(sequence,Metal=Metal,NbindingMetal=NbindingMetal,Maxcharge=Maxcharge,
                                                    Adduct=Adduct
    )
  }
  ESI.Theoretical_res<-ESI.Theoretical(sequence,Metal=Metal,NbindingMetal=NbindingMetal,Maxcharge=Maxcharge,
                                       Adduct=Adduct
  )
  ESI.Theoretical_res
  TargetESI_results<-TargetESI(parent=parent,ESI.Theoretical_res,mzw = mzw ,Adduct=Adduct)
  colnames(TargetESI_results)<-c("Monoisotopic mass","Mass error","Stoichiometry","Charge","Adduct")
  write.csv2(TargetESI_results,paste("CID_precursor_",parent,".csv",sep=""))

  Metals<-TargetESI_results[,3]

  Nmetal<-as.numeric(regmatches(Metals, gregexpr("[[:digit:]]+", Metals))[[1]][1])


  theoryMultiple<-MS.MS.Multiple(sequence,native=TRUE,Metal=Metal,Nmetal=NbindingMetal,Adduct=Adduct)


  if (peak==TRUE){
    exp<-experimental
    exp <-exp[which(exp[, 1] > mz1& exp[, 1] < mz2) , ]
    peakInfo <- peakDetectionCWT(exp[,2],SNR.Th=SNR.Th,nearbyPeak=TRUE)
    majorPeakInfo = peakInfo$majorPeakInfo
    peakIndex <- majorPeakInfo$peakIndex
    betterPeakInfo <- tuneInPeakInfo(exp[,2], majorPeakInfo)


    mz<- exp[,1][betterPeakInfo$peakIndex]

    DetectedSmoothed<-cbind(mz,betterPeakInfo$peakValue)

  MassSpecWavelet::plotPeak(exp[,2],betterPeakInfo$peakIndex
    )
    readline(prompt="Peak detection: Press [enter] to accept it and continue")

    #plot(DetectedSmoothed,type="h")
    experimental<-DetectedSmoothed
  }

  Deconv_CysMpro<-Deconvolved_ESI_High(exp=experimental,maxMass=maxMass,SNR.Th=SNR.Th,max.charge=Maxcharge,tol=tol)

  write.csv2(Deconv_CysMpro[[1]],paste("Deconv_CysMpro_all",parent,".csv",sep=""))

  write.csv2(Deconv_CysMpro[[2]],paste("Deconv_CysMpro_",parent,".csv",sep=""))

  final<-target.MSMS.Multiple(theoryMultiple,Deconv_CysMpro,t=t,b=b,native=TRUE,plot= plot)

  #re<-final[order(final$expt_mz),]

  return(list(final,Deconv_CysMpro))
}
