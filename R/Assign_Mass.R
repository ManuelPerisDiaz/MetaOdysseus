#' Annotate ESI-MS spectra
#'
#' Annotate ESI-MS spectra

#' @param sequence Protein sequence.
#' @param Metal Metal ion to include in the calculation of the MS library.
#' @param NbindingMetal Number of expected metal ions binding to the protein. Default to 7.
#' @param Scored Scored ESI-MS spectra
#' @param mono If TRUE, monoisotopic mass is used.
#' @param mzw Mass error window
#' @param write If TRUE, reports the results in csv file
#' @param dir Directory where to save the results
#' @param Maxcharge Maximum charge expected for adduct considered. Default to 6.
#' @param Adduct Adduct for library generation.
#' @return Annotate ESI-MS spectra
#' @export

Assign_Mass<- function(Scored,sequence, Metal, mono=FALSE, NbindingMetal = 7,Maxcharge=6,Adduct=NULL, mzw = 1, write=FALSE,dir=getwd()){


  peptideVector <- strsplit(sequence, split = "")[[1]]
  resultsVector <- c(C = 0, H = 0, N = 0, O = 0, S = 0)
  for (i in 1:length(peptideVector)) {
    resultsVector <- FindElement(peptideVector[i]) + resultsVector
  }

  resultsVector <- resultsVector + c(C = 0, H = 2, N = 0, O = 1, S = 0)

  Zn=63.928047
  Cd=113.903358
  Hg=201.970596
  Cu=62.929600
  Pb=207.976593
  NH4=18.033823
  Na=23
  K=39
  Fe=55.934940



  if (is.null(Metal)){
    element<-NULL
  }else{
    if (Metal== "Zn"){
      element<-Zn}
    if (Metal=="Cd"){
      element<-Cd}
    if (Metal=="Hg"){
      element<-Hg}
    if (Metal=="Cu"){
      element<-Cu}
    if (Metal=="Pb"){
      element<-Pb}
    if (Metal=="Fe"){
      element<-Fe}

  }



  massH<-c()
  nrow = as.numeric(NbindingMetal)+1
  massH<-matrix(nrow = as.numeric(NbindingMetal)+1)





    for (i in 1:nrow(massH)){

      if(mono==TRUE){
      massH[i,] <- OrgMassSpecR::MonoisotopicMass(as.list(c(resultsVector,x=i-1)),isotopes=list(x=element),charge=0)
      }else{
        massH[i,] <-OrgMassSpecR::MolecularWeight(as.list(c(resultsVector,x=i-1)),amu=list(x=element))
      }



      MASSlist<- do.call(rbind,lapply(massH,matrix,ncol=1,byrow=TRUE))

  }
  rownames(massH)<-paste0("Metal",1:nrow-1,sep="")


  exp<-unlist(Scored$Molecular_weight)
  Uscored<-unlist(Scored$Uscore)
  Dscored<-unlist( Scored$DScore)
  Fscored<-unlist( Scored$Fscore)
  CCScored<-unlist( Scored$CCScore)
  MScored<-unlist( Scored$MScore)


  databaseH<-as.data.frame(massH)
  names<-row.names(massH)

  f<-mzw
  ions1bb<-data.frame()
  msbc<-data.frame()
  b<-data.frame()
  msb<-data.frame()
  ms <- data.frame()


  dat <- as.data.frame(matrix(ncol=8, nrow=0))


  for (i in 1:nrow(databaseH)){

    fullmsset <- databaseH[i, 1]
    options(digits=6)
    fullms2 <-  as.numeric(fullmsset[!is.na(fullmsset)])
    name<- names[i]



     msb <- subset(exp, exp >  fullms2 - f & exp <  fullms2 + f )
    Indexa<-match(msb,exp)

    Uscored_<-Uscored[Indexa]
    Dscored_<-Dscored[Indexa]
    Fscored_<-Fscored[Indexa]
    CCScored_<- CCScored[Indexa]
    MScored_<- MScored[Indexa]

    if (length(msb) > 0){
        ppmdif<-((msb[1])- fullms2)

         dat[i,1]<-name
        dat[i,2]<-msb
        dat[i,3]<- ppmdif
        dat[i,4]<-Dscored_
          dat[i,5]<-Uscored_
          dat[i,6]<- Fscored_
          dat[i,7]<-   CCScored_
          dat[i,8]<- MScored_
        }


  }
dat<-  na.omit(dat)
  colnames(dat)<-c("Protein","Mass","Error (Da)","Dscore","Uscore","Fscore","CCScore","MScore")

  dat <- dat[dat$Uscore >= 0, ]

      if(write==TRUE){
        options(digits = 2)
        write.csv2(dat,paste0("ESIannotated",".csv") )
      }

      return(dat)

    }



