

CountAAPMF<-function(seq,maldiTOF){

ff<-c(maldiTOF$Seq)
mt3clv<-Digest(seq, missed=2,enzyme = "trypsin",IAA=FALSE)


Results<-c()
matched2<-c()
for (i in 1:length(ff)){
  matched<-c()

 i<- match(ff[i],mt3clv[,1] )
  matched <-mt3clv[i,]
matched2<-rbind(matched2,matched)
}


Results<-cbind(maldiTOF,matched2[,2:4])
Resultssubseted<-Results[!duplicated(Results[6]),]
total<-as.numeric()

for (i in 1:length(Resultssubseted$Seq)){
  #for (i in 1:length(unique(Results$Seq))){
#for (i in 1:length((Results$Seq))){
  peptide_vector <- strsplit(Resultssubseted$Seq, split = "")[[i]]
  #peptide_vector <- strsplit(unique(Results$Seq), split = "")[[i]]
# peptide_vector <- strsplit((Results$Seq), split = "")[[i]]

peptide_length <- length(peptide_vector)
total<-sum(as.numeric(peptide_length),total)
}



########### for the peptides identifiy, make the interval numeric, count, sum the total number,
########### and optional is to check if some AA has been sumed several times because overlapping
########## of sequences, check how many overlappings and then rest to the sum previously done

TotalOverlappingsb<-c()
sed22<-c()
RR<-matrix(,5000,5)
for (i in 1:nrow(Resultssubseted)){
   sed=NULL
   RR<-matrix(,5000,5)
   sed<- as.numeric(Resultssubseted[i,7:8])
   Max<- max(nrow(Resultssubseted))


      sed2=seq2=seq=NULL

          for (j in 1:Max){
            seq<-as.numeric(Resultssubseted[j,7:8])
            sed2<-seq( sed[1], sed[2],by=1)
            seq2<-seq(seq[1],seq[2],by=1)

            overl<-as.numeric(which(sed2%in%seq2))
            overllop<-length(overl)
            RR[j,1]<-seq[1]
            RR[j,2]<-seq[2]
            RR[j,3]<-sed[1]
            RR[j,4]<-sed[2]
            RR[j,5]<-overllop

               }
      FinalR<-RR[!is.na(RR[,5]),]
      TotalOverlappings<-c()
      if (NCOL(FinalR)>1){
      TotalOverlappingsC<-(sum(FinalR[,5]))}else{ TotalOverlappings<-((FinalR[5]))


      }

      NoOverlaping<-subset(FinalR,FinalR[,5]==0)


      TotalOverlappings<-TotalOverlappingsC-length(sed2)
      TotalOverlappingsb<-cbind(TotalOverlappingsb,TotalOverlappings)


}

sum(TotalOverlappingsb)

FinalR<-RR[!is.na(RR[,5]),]

TotalOverlappings<-c()
TotalOverlappings<-(sum(FinalR[,5]))


AaCountedRested<-total-TotalOverlappings

return(AaCountedRested)
}
