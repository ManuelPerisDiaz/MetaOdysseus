#' Sequence coveraged from Bottom-up MALDI-MS experiment.
#'
#' Calculates the protein sequence coverage provided by bottom-up MALDI-MS experiments
#'
#' @param seq Protein sequence.
#' @param maldiTOF Output obtained from \code{\link{annotationMALDI_Bot}}.
#' @param missed Missed cleavage sites produced by proteolytic enzyme. Default to 2.
#' @param enzym Proteolytic enzyme used in the experiment. Default to Trypsin.
#' @examples
#' \dontrun{
#' MT2a<-("MDPNCSCAAGDSCTCAGSCKCKECKCTSCKKSCCSCCPVGCAKCAQGCICKGASDKCSCCA")
#'TheoreticalMass<-Protein.to.Peptide(MT2a)
#'ProcessMALDIMS<-ProcessMALDIMS(file,t1=500,t2=3500,plotProcess=TRUE,plotDetection=TRUE,plotOriginal = TRUE,comb = "comb3")
#'maldiTOF<-annotationMALDI_Bot(ProcessMALDIMS,TheoreticalMass,mzw = 1,plot=TRUE,plotIndi=TRUE,orig=TRUE,write=TRUE,In=10)
#' maldiTOF.Amp(MT2a,maldiTOF)
#' CountAAPMFF<-SequenceCounting(MT2a,maldiTOF)
#' }
#' @return Protein sequence coverage for bottom-up MALDI-MS experiment.
#' @export

SequenceCounting<-function(seq,maldiTOF,missed,enzym = "trypsin"){

  Digest  = getFromNamespace("Digest", "OrgMassSpecR")

  peptide_vector <- strsplit(seq, split = "")[[1]]
  peptide_length <- length(peptide_vector)
   ff<-c(maldiTOF$Seq)
   mt3clv<-Digest(seq=seq, missed= missed,enzym=enzym,IAA=FALSE)


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
  Resultssubseted<-  Resultssubseted[complete.cases(Resultssubseted), ]


 # if (nrow(  Resultssubseted)>0){

RR<-matrix(,5000, peptide_length+10 )
for (i in 1:nrow(Resultssubseted)){

   sed2<-seq( Resultssubseted[i,7], Resultssubseted[i,8],by=1)
   from<-min(sed2)
     to<-max(sed2)
   RR[i,from:to]<-sed2

}
tm1<-RR[,-which(apply(RR,2,function(x)all(is.na(x))))]
tm2<-tm1[-which(apply(RR,1,function(x)all(is.na(x)))),]

vector<-matrix(,1, peptide_length )


SeqCoverage <-(length(tm2)/ peptide_length )*100

  #}
#else{

 #   SeqCoverage <-0}

return(SeqCoverage)

}


