#' Annotate top-down MS/MS experiments from multiple Cys-protein.
#'
#' @param theoryMultiple Output obtain from \code{\link{MS.MS.Multiple}}.
#' @param exp Experimental MS/MS file.
#' @param t Mass-to-charge window (Da) for the annotation. Default to 0.4.
#' @param b Threshold for relative intensity. Default to 5.
#' @param supress Include MS/MS ions in the output plot. Default to FALSE.
#' @param xlim Mass range considered for the annotation.
#' @param write Output the results
#' @param plot plot the MS/MS spectra
#' @examples
#' @return Peaks annotated in the top-down MS/MS spectra.
#' @export

target.MSMS.Multiple_B<-function(theoryMultiple,exp,t=0.4,b=5,write=FALSE,plot=FALSE,supress = FALSE,xlim = NULL){

  if (is.null(xlim)){
  t1<- min(exp[,1])
  t2<-max(exp[,1])

  xlim = c(t1, t2)
  }


  final<-c()
  for (i in 1:length(theoryMultiple)){
    results<-c()

    if(plot==TRUE){
    tiff(paste0("theoryMultiple",i,".tiff",sep=""),width = 9, height = 6, units = 'in',res=600)
    }
       results<-targetMSMS_bottom(Deconv_CysMpro=exp,theory=theoryMultiple[[i]],plot=plot,mzw=t,b=b,

                           xlim =xlim,supress = supress)



       if(plot==TRUE){
        dev.off()}
  final[i]<-list(results)


if(write==TRUE){
  write.csv2(final[i],paste0("Sequence",i,".csv",sep=""))}

   }

  return(final)
}


