#' Annotate a set of top-down MS/MS spectra
#'
#' @param t Mass error window.
#' @param mz1 left mass
#' @param mz2 right mass
#' @param Maxcharge maximum charge state for deconvolution
#' @param tol mass tolerance
#' @param b Filter by intensity
#' @param SNR.Th SNR threshold
#' @param peak if TRUE peak picking is performed
#' @param plot if TRUE, plot is produced
#' @param maxMass maximum mass for deconvolution
#' @param mzw mass error window to annotate the parent ion
#' @param NbindingMetal Number of metal ions
#' @param Metal Metal
#' @param label Label the results
#' @param mzML if TRUE, the experimental file is in mzML file format
#' @param txt if TRUE, the experimental file is in txt file format
#' @param xy if TRUE, the experimental file is in xy file format
#' @param sequence Protein sequence
#' @param dir Where to save results
#' @examples
#' @return Annotate a set of top-down MS/MS spectra
#' @export

Topdown_batch<-function(mzML=TRUE,xy=FALSE,mz1=500,mz2=3000,tol=0.15,

                        txt=FALSE,Metal=c("Zn"),SNR.Th=2,peak=FALSE,
                        plot=FALSE,dir=getwd(),sequence, b,t,maxMass,mzw,
                        NbindingMetal,Maxcharge,label="CID"){


  if (mzML==TRUE){
    folder<-dir()[grep(".mzML",dir())]

  }
  if (txt==TRUE){
    folder<-dir()[grep(".txt",dir())]
  }
  if (xy==TRUE){
    folder<-dir()[grep(".xy",dir())]
  }

  matches_f <- vector("list")
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

    parentb<-parent

    result<- Topdown(experimental=exp,parent=parent,SNR.Th=SNR.Th,peak=peak, mz1=mz1, mz2=mz2, tol=tol,sequence=sequence,t=t, b=b,plot=plot,maxMass=maxMass ,mzw=mzw, Metal=Metal, NbindingMetal=NbindingMetal,Maxcharge=Maxcharge,Adduct=NULL)
    write.csv2(result,paste(parent[1],label,".csv",sep=""))

    new.folder <-paste(parent[1],label,sep="")
    file.copy(folder[z], new.folder)
    file.copy(paste(parentb[1],label,".csv",sep=""), new.folder)
    file.copy(paste("CID_precursor_",parentb[1],".csv",sep=""), new.folder)
    file.copy(paste("Deconv_CysMpro_",parentb[1],".csv",sep=""), new.folder)
    file.copy(paste("Deconv_CysMpro_all",parentb[1],".csv",sep=""), new.folder)


    matches_f[[z]]<-result

  }
  return( matches_f)
}
