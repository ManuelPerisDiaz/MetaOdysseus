
ESI.theoretical.alkylators<-function (sequence, Metal=c("Zn"), NbindingMetal,Maxcharge,Adduct=c("NH4"),Alkylator=c("IAM"))
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


  M<- c(C = 0, H = 0, N = 0, O = 0, S = 0)
  M2<-c(C = 0, H = 0, N = 0, O = 0, S = 0)

  for (i in 0:sum(as.numeric(peptideVector=="C"))){
    # formula for NEM
    if (Alkylator=="NEM"){
    elementx<-c(C=6,H=7,N=1,O=2,S=0)}
    if (Alkylator=="IAM"){
      elementx<-c(C=2,H=3,N=1,O=1,S=0)}

    M<- resultsVector+elementx*i
    M2[i+1]<-list(M)
  }


   mod= sum(as.numeric(peptideVector=="C"))

   nrow = as.numeric(NbindingMetal)+1


   withmod<-matrix(nrow = (length(M2)-1), ncol = as.numeric(Maxcharge))

   massH<-matrix(nrow =  nrow, ncol = as.numeric(Maxcharge))

 for (b in 1:length(M2)){
   for (i in 1:(as.numeric(NbindingMetal)+1)){
    for (a in 1:Maxcharge){

       massH[i,a] <- OrgMassSpecR::MonoisotopicMass(as.list(c(M2[[b]],x=i-1)),isotopes=list(x=element),charge=a)

       rownames(massH)<-paste0("Metal",1:nrow-1,sep="_","Mod",b-1)

      colnames(massH)<-paste0("Charge",1:Maxcharge,sep="")



    }
   }

   withmod<- rbind(withmod,massH)

   }
   withmod <- na.omit(withmod)
   MASSlist<- do.call(rbind,lapply( withmod,matrix,ncol=1,byrow=TRUE))

   massH <- withmod

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
