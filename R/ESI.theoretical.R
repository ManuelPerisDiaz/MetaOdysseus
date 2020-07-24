#'  ESI-MS library construction for metal-binding protein
#'
#' Generates a theoretical ESI-MS library for a particular metal binding protein selecting the desired metal ion, metal ion stoichiometry and adducts.
#'
#' @param sequence Protein sequence.
#' @param Metal Metal ion to include in the calculation of the MS library.
#' @param NbindingMetal Number of expected metal ions binding to the protein. Default to 7.
#' @param Maxcharge Maximum charge expected for adduct considered. Default to 6.
#' @param Adduct Adduct for library generation.
#' @return ESI-MS library for metal-binding protein.
#' @examples
#' \dontrun{
#' MT2a<-("MDPNCSCAAGDSCTCAGSCKCKECKCTSCKKSCCSCCPVGCAKCAQGCICKGASDKCSCCA")
#'Theoretical.esi<-ESI.Theoretical(MT2a, Metal=c("Zn"), NbindingMetal = 7,Maxcharge=6,Adduct=c("Na"))
#' }
#' @export

ESI.Theoretical<-function (sequence, Metal=c("Zn"), NbindingMetal,Maxcharge,Adduct=c("NH4"))
{
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

charges<-c()
  nrow = as.numeric(NbindingMetal)+1
  massH<-matrix(nrow = as.numeric(NbindingMetal)+1, ncol = as.numeric(Maxcharge))
  for (i in 1:nrow(massH)){
    for (a in 1:Maxcharge){

      massH[i,a] <- OrgMassSpecR::MonoisotopicMass(as.list(c(resultsVector,x=i-1)),isotopes=list(x=element),charge=a)
      rownames(massH)<-paste0("Metal",1:nrow-1,sep="")
      colnames(massH)<-paste0("Charge",1:Maxcharge,sep="")
      MASSlist<- do.call(rbind,lapply(massH,matrix,ncol=1,byrow=TRUE))
    }
  }

  if (is.null(Adduct)){

    return(list(massH))}else{



  massNH4<-c()
if (Adduct== "NH4"){
  massNH4<-massH+NH4
  return(list(massH,massNH4))

}

  massNa<-c()
  if (Adduct== "Na"){
    massNa<-massH+Na
    return(list(massH,massNa))

  }

  massK<-c()
  if (Adduct== "K"){
    massK<-massH+K
    return(list(massH,massK))

  }


    }
}
