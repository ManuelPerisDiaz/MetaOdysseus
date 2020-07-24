
Protein.To.Peptide_dual<-function (sequence,missed=2,enzym = "trypsin")
{{


  mt3clv<-Digest(sequence, enzym ,missed,IAA=FALSE)

  result<-c()
  result<-vector("list",length(mt3clv[,1]))



  for (i in 1:length(mt3clv[,1])){
    seq<-mt3clv[i,1]
    name<-seq
    num<-i
    peptideVector <- strsplit(seq, split = "")[[1]]
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

    NCys<-sum(as.numeric(peptideVector=="C"))
    M2<-c()
    M<-c()
    M4<-c()
    M5<-c()
    M6<-c()
    for (i in 1:NCys){
      element<-c(C=2,H=3,N=1,O=1,S=0)
      M<- resultsVector+element*i
      M2[i]<-list(M)}

      for (i in 1:length(M2)){
        M2s<-M2[[i]]
        nNEM<- NCys-i
        NEM<-c(C=6,H=6,N=1,O=2,S=0)


          for (i in 0:nNEM){

            M3<-  M2s+NEM*i
             M4[i+1]<-list(M3)

             }



        M6<-list(M5)

        }


    }

## we need to add when NEM=0 as well, only consider now with NEM.


    mass<-c()
    massb<-c()
    for (i in 1:length(M2)){
      mass[i]<- MonoisotopicMass(M2[[i]],charge=1)
    }
    i<-num
    hist(mass)
    NM<-MonoisotopicMass(resultsVector)
    Modifications=seq(0,NCys,1)
    mass<-c(NM,mass)
    massb<-cbind(Modifications=Modifications,mass)

    result[[i]]<-massb
    names(result)[[i]]<-paste(name)

  }

  return(  result)

}

