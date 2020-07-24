#' Library generation of chemical protein Cys-labelled residues.
#'
#' Generates theoretical list of Cys-labelled residues with a particular chemical labelling reagent.
#'
#' @param sequence Protein sequence.
#' @param IAA Selection of chemical labelling reagent used: Iodoacetic acid.
#' @param IAM Selection of chemical labelling reagent used: Iodoacetamide.
#' @param ET Selection of chemical labelling reagent used: Ethyl iodoacetate.
#' @param NEM Selection of chemical labelling reagent used: NEM.
#' @param LABEL Chemical labelling. Default to TRUE.
#' @param monoisotopic Monoisotopic mass considered for library construction. Default to FALSE.
#' @examples
#' \dontrun{
#' MT2a<-("MDPNCSCAAGDSCTCAGSCKCKECKCTSCKKSCCSCCPVGCAKCAQGCICKGASDKCSCCA")
#'TheoreticalMass<-Protein.to.List(MT2a, ET= TRUE)
#' @return MS library generated.
#' @export


Protein.to.List<-function (sequence, IAA = FALSE,IAM=FALSE, ET=FALSE, NEM=FALSE,LABEL=TRUE,monoisotopic=FALSE)
{
  peptideVector <- strsplit(sequence, split = "")[[1]]

    FindElement <- function(residue) {
      if (residue == "A")
        element <- c(C = 3, H = 5, N = 1, O = 1, S = 0)
      if (residue == "R")
        element <- c(C = 6, H = 12, N = 4, O = 1, S = 0)
      if (residue == "N")
        element <- c(C = 4, H = 6, N = 2, O = 2, S = 0)
      if (residue == "D")
        element <- c(C = 4, H = 5, N = 1, O = 3, S = 0)
      if (residue == "E")
        element <- c(C = 5, H = 7, N = 1, O = 3, S = 0)
      if (residue == "Q")
        element <- c(C = 5, H = 8, N = 2, O = 2, S = 0)
      if (residue == "G")
        element <- c(C = 2, H = 3, N = 1, O = 1, S = 0)
      if (residue == "H")
        element <- c(C = 6, H = 7, N = 3, O = 1, S = 0)
      if (residue == "I")
        element <- c(C = 6, H = 11, N = 1, O = 1, S = 0)
      if (residue == "L")
        element <- c(C = 6, H = 11, N = 1, O = 1, S = 0)
      if (residue == "K")
        element <- c(C = 6, H = 12, N = 2, O = 1, S = 0)
      if (residue == "M")
        element <- c(C = 5, H = 9, N = 1, O = 1, S = 1)
      if (residue == "F")
        element <- c(C = 9, H = 9, N = 1, O = 1, S = 0)
      if (residue == "P")
        element <- c(C = 5, H = 7, N = 1, O = 1, S = 0)
      if (residue == "S")
        element <- c(C = 3, H = 5, N = 1, O = 2, S = 0)
      if (residue == "T")
        element <- c(C = 4, H = 7, N = 1, O = 2, S = 0)
      if (residue == "W")
        element <- c(C = 11, H = 10, N = 2, O = 1, S = 0)
      if (residue == "Y")
        element <- c(C = 9, H = 9, N = 1, O = 2, S = 0)
      if (residue == "V")
        element <- c(C = 5, H = 9, N = 1, O = 1, S = 0)
      if (residue == "C")
        element <- c(C = 3, H = 5, N = 1, O = 1, S = 1)
      return(element)
    }
    resultsVector <- c(C = 0, H = 0, N = 0, O = 0, S = 0)

    for (i in 1:length(peptideVector)) {
      resultsVector <- FindElement(peptideVector[i]) +
        resultsVector
    }
    resultsVector <- resultsVector + c(C = 0, H = 2, N = 0,
                                       O = 1, S = 0)

    if (IAM== "TRUE"){
    M<- c(C = 0, H = 0, N = 0, O = 0, S = 0)
    M2<-c(C = 0, H = 0, N = 0, O = 0, S = 0)

    for (i in 1:sum(as.numeric(peptideVector=="C"))){
    element<-c(C=2,H=3,N=1,O=1,S=0)
    M<- resultsVector+element*i
    M2[i]<-list(M)
    }

    mass<-c()

    if (monoisotopic=="TRUE"){
      for (i in 1:length(M2)){
        mass[i]<- MonoisotopicMass(M2[[i]],charge=1)
      }
      hist(mass)

      NM<-MonoisotopicMass(resultsVector,charge=1)
    }else{
      for (i in 1:length(M2)){
        mass[i]<- AverageMass(M2[[i]],charge=1)
      }
      hist(mass)

      NM<-AverageMass(resultsVector,charge=1)

    }

    Modifications=seq(0,sum(as.numeric(peptideVector=="C")),1)

    mass<-c(NM,mass)
    mass<-cbind(Modifications=Modifications,mass)
    return(( mass))


    }

    if (NEM== "TRUE"){
      M<- c(C = 0, H = 0, N = 0, O = 0, S = 0)
      M2<-c(C = 0, H = 0, N = 0, O = 0, S = 0)

      for (i in 1:sum(as.numeric(peptideVector=="C"))){
        element<-c(C=6,H=7,N=1,O=2,S=0)
        M<- resultsVector+element*i
        M2[i]<-list(M)
      }

      mass<-c()

      if (monoisotopic=="TRUE"){
      for (i in 1:length(M2)){
        mass[i]<- MonoisotopicMass(M2[[i]],charge=1)
      }
      hist(mass)

      NM<-MonoisotopicMass(resultsVector,charge=1)
      }else{
        for (i in 1:length(M2)){
          mass[i]<- AverageMass(M2[[i]],charge=1)
        }
        hist(mass)

        NM<-AverageMass(resultsVector,charge=1)
      }

      Modifications=seq(0,sum(as.numeric(peptideVector=="C")),1)

      mass<-c(NM,mass)
      mass<-cbind(Modifications=Modifications,mass)
      return(( mass))
    }

    if (IAA== "TRUE"){
      M<- c(C = 0, H = 0, N = 0, O = 0, S = 0)
      M2<-c(C = 0, H = 0, N = 0, O = 0, S = 0)

      for (i in 1:sum(as.numeric(peptideVector=="C"))){
        element<-c(C=2,H=2,N=0,O=2,S=0)
        M<- resultsVector+element*i
        M2[i]<-list(M)
      }

      mass<-c()
      if (monoisotopic=="TRUE"){
        for (i in 1:length(M2)){
          mass[i]<- MonoisotopicMass(M2[[i]],charge=1)
        }
        hist(mass)

        NM<-MonoisotopicMass(resultsVector,charge=1)
      }else{
        for (i in 1:length(M2)){
          mass[i]<- AverageMass(M2[[i]],charge=1)
        }
        hist(mass)

        NM<-AverageMass(resultsVector,charge=1)

      }

      Modifications=seq(0,sum(as.numeric(peptideVector=="C")),1)

      mass<-c(NM,mass)
      mass<-cbind(Modifications=Modifications,mass)
      return(( mass))
    }


    if (ET== "TRUE"){
      M<- c(C = 0, H = 0, N = 0, O = 0, S = 0)
      M2<-c(C = 0, H = 0, N = 0, O = 0, S = 0)

      for (i in 1:sum(as.numeric(peptideVector=="C"))){
        element<-c(C=4,H=6,N=0,O=2,S=0)
        M<- resultsVector+element*i
        M2[i]<-list(M)
      }

      mass<-c()
      if (monoisotopic=="TRUE"){
        for (i in 1:length(M2)){
          mass[i]<- MonoisotopicMass(M2[[i]],charge=1)
        }
        hist(mass)

        NM<-MonoisotopicMass(resultsVector,charge=1)
      }else{
        for (i in 1:length(M2)){
          mass[i]<- AverageMass(M2[[i]],charge=1)
        }
        hist(mass)

        NM<-AverageMass(resultsVector,charge=1)

      }

      Modifications=seq(0,sum(as.numeric(peptideVector=="C")),1)

      mass<-c(NM,mass)
      mass<-cbind(Modifications=Modifications,mass)
      return(( mass))
    }

    if (LABEL== "FALSE"){
    return(as.list(resultsVector))
    }
    }





