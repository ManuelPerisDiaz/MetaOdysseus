
SNR_loop<-function(chargemax=50,maxMass=7000,W=1, SNRuser=20,ordered.mz.snr,hmz){

SNRloop=0
SNRloopr=0
chargemax=chargemax
maxMass<-maxMass
W=W
SNRuser=SNRuser
SNRloop=0
SNRloopr=0

dat <- as.data.frame(matrix(ncol=3, nrow=0))
M3<- as.data.frame(matrix(ncol=3, nrow=0))

b <- c()

vectorCharge<-c(1:chargemax)
for(i in 1:length(vectorCharge)){
  charge<-vectorCharge[i]
  m=charge*hmz
  remainingcharge<-vectorCharge[-i]
  calc.mz<-m/remainingcharge
  semi<-cbind(charge=remainingcharge,mz=calc.mz)
  full<-cbind(charge=charge,mz=hmz)
  fullx<-rbind(full,semi)


  totalscore=0
  score=0

  charge1<-fullx[which(fullx[, 1] ==1 ),]
  Mz<-as.numeric(charge1[2])

  if( Mz<=maxMass){
  for(i in 1:length(fullx[,1])){
    calc.mz<-as.numeric(fullx[i,2])
    exp.mz <- ordered.mz.snr[which(ordered.mz.snr[, 1] > calc.mz-W & ordered.mz.snr[, 1] < calc.mz+W),
                             ]

    while(nrow(exp.mz)>1 ){
      W<-W-0.01
      exp.mz <- ordered.mz.snr[which(ordered.mz.snr[, 1] > calc.mz-W & ordered.mz.snr[, 1] < calc.mz+W),
                               ]
    }

    if(nrow(exp.mz)==1){

      if( exp.mz[,2]> SNRuser){

        score<-log(exp.mz[,2])
        dat[i,1]<-exp.mz[,1]
        dat[i,2]<-score
        dat[i,3]<-exp.mz[,2]
        dat<-dat[!duplicated(dat),]
        dat<-dat[complete.cases(dat),]
        totalscore<-totalscore+score

      }else{
        score<-log(exp.mz[,2])
        M3[i,1]<-exp.mz[,1]
        M3[i,2]<-score
        M3[i,3]<-exp.mz[,2]
        M3<-M3[!duplicated(M3),]
        M3<-M3[complete.cases(M3),]
      }

    }


  }
  if (nrow(M3)>0){
    SNRloop<-M3$V3
    SNRloopr<-c(SNRloopr, SNRloop)

  }


  SNRthreshold<-max(SNRloopr)
  if(SNRthreshold==0){
    SNRthreshold<-0.01

  }


  }

  }
return(SNRthreshold)
}
