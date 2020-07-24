#' Annotate Bottom-up MALDI-MS experiments
#'
#' Annotate Bottom-up MALDI-MS experiments
#'
#' @param ProcessedMALDIMS Output obtain from \code{\link{ProcessMALDIMS}}.
#' @param TheoreticalMass Output obtained from \code{\link{Protein.to.Peptide}}.
#' @param mzw Mass-to-charge window (Da) for the annotation. Default to 1.
#' @param write Write an output table with results. Default to  FALSE.
#' @param plot Plot the peptides found. Default to TRUE.
#' @param orig Plot raw MS spectra. Default to TRUE.
#' @param lint Left interval for plotting purposes. Default to 3.
#' @param rint Right interval for plotting purposes. Default to 3.
#' @param FilterInt Control filtering by relative intensity. Default to TRUE.
#' @param In value for filtering by relative intensity. Default to 10.
#' @param plotIndi Plot specific peptides matched. Default to TRUE.
#' @param dir Directory where to save plots. Default to C:.
#' @param save.name Name to save output.
#' @examples
#' @return Peaks annotated in the MALDI-MS spectra.
#' @export


annotationMALDI_Bot<- function (ProcessedMALDIMS,TheoreticalMass, mzw = 1, write=FALSE,plot=TRUE,orig=TRUE,
                                lint=3,rint=3,FilterInt=TRUE,In=10,plotIndi=TRUE,dir=getwd(),save.name

                                ){
  exp=NULL
  or=NULL
  exp<-as.data.frame(ProcessedMALDIMS$IN.Mz)
  or<-as.data.frame(ProcessedMALDIMS$OR)
  y<-exp[,2]
  yx = (y-min(y))/(max(y)-min(y))*100
  exp<-cbind(exp,ppmdif="ppmdif")
  exp<-cbind(exp,N.Modifications="N.Modifications")
  exp<-cbind(exp,yx)
  msbcc<-data.frame()

  for (i in 1:length(TheoreticalMass)){


    nameTry<-names(TheoreticalMass[i])

    database<-TheoreticalMass[[i]]



  f<-mzw
  ions1bb<-data.frame()
  msbc<-data.frame()
  msb<-data.frame()
  b<-data.frame()
  ms <- data.frame()
  for (i in 1:nrow(database)){
    from<-2
    to<-as.numeric(length(database[i,]))
    fullmsset <- database[i, from:to]
    options(digits=6)
    fullms2 <-  as.numeric(fullmsset[!is.na(fullmsset)])
    name<- database[i, 1]


    for (i in 1:length(fullms2)){
      msb <- subset(exp, exp[, 1] >  fullms2[i] - f & exp[, 1] <  fullms2[i] + f )
      if (nrow(msb) > 0){
        ppmdif<-((msb[1])- fullms2[i])
        msb[,3]<- ppmdif

        msb[,4]<-name
        msb[,6]<-nameTry

      }
      msbc<-rbind(msb,msbc)
      msbc<-as.data.frame(msbc)

    }

  }

  if(nrow(msbc)>0){
  colnames(msbc) <- c("m/z","Maximum Intensity","Error(Da)","N.Modifications","Rel.Intensity","Seq")

  msbcc<-rbind(msbcc,msbc)
  msbcc<-as.data.frame(msbcc)

  colnames(msbcc) <- c("m/z","Maximum Intensity","Error(Da)","N.Modifications","Rel.Intensity","Seq")




  if(FilterInt==TRUE){

    msbcc<-subset(msbcc, msbcc[,5]>In)

  }

   if(plot==TRUE & nrow(msbc)>0){

    t1<-min(msbc[,1])-lint
    t2<-max(msbc[,1])+rint
    plotPeakBot(ProcessedMALDIMS,t1=t1,t2=t2,orig=orig,nameTry=nameTry,dir=dir)}

  }

  if(plotIndi==TRUE){

    for (i in 1:nrow(msbcc) ){
    mass<-msbcc[i,]$`m/z`
   subSe<- subset(ProcessedMALDIMS$OR, ProcessedMALDIMS$OR$peak...1. >  msbcc[i,]$`m/z` - 1 &
             ProcessedMALDIMS$OR$peak...1. <  msbcc[i,]$`m/z` + 3 )
   plotPMFPeak(subSe,orig = FALSE,norm=TRUE,QC1=FALSE,QC2=FALSE,mass=mass,dir=dir)


  }}


  if(write==TRUE ){
    write.csv2(msbcc,paste0("MALDIannotated",save.name,".csv"))}
  }

  if(nrow(msbcc)>0){
  return(msbcc)}else{print("No matched peaks")}
}






