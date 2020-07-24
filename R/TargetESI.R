#' Annotate ESI-MS spectra
#'
#' @param parent Parent mass
#' @param Theoretical.esi Theoretical library
#' @param Adduct Adduct
#' @param mzw Mass-to-charge window (Da) for the annotation. Default to 1.
#' @examples
#' @return Peaks annotated in the ESI-MS spectra.
#' @export

TargetESI<- function (parent,Theoretical.esi,Adduct, mzw = 1){
  exp=NULL
  or=NULL
  exp<-as.data.frame(parent)
  or<-as.data.frame(Theoretical.esi)

  if (!is.null(Adduct)){
    databaseH<-data.frame(Theoretical.esi[[2]])
  }else{

    Adduct=c("H+")
  databaseH<-data.frame(Theoretical.esi[[1]])
  }
  f<-mzw
  ions1bb<-data.frame()
  msbc<-data.frame()
  b<-data.frame()

  ms <- data.frame()

  msb<-matrix(ncol= 5)


  for (i in 1:nrow(databaseH)){
    from<-1
    to<-as.numeric(length(databaseH[i,]))
    fullmsset <- databaseH[i, from:to]
    options(digits=6)
    fullms2 <-  as.numeric(fullmsset[!is.na(fullmsset)])
    name<- row.names(databaseH[i,])

    for (i in 1:length(fullms2)){
      ms <- exp[which(exp[, 1] >  fullms2[i] - f & exp[, 1] <  fullms2[i] + f ),]
      Chargename<- colnames(databaseH[i])
      if (length(ms) > 0){


         ppmdif<-((ms[1])- fullms2[i])

         msb[,1] <-  ms
         msb[,2]<- ppmdif

        msb[,3]<-name
        msb[,4]<-Chargename
        msb[,5]<- Adduct
      }
       }
    }

return(msb)
}


