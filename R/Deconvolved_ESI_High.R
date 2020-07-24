#' Pre-process and deconvolved ESI-MS spectra with Zscore algorithm
#'
#' @param exp Experimental MS file.
#' @param mz1 Select lower m/z limit.
#' @param SNR.Th Signal-to-noise ratio. Default to 3.
#' @param max.charge Maximum charge state
#' @param maxMass Maximum mass deconvoluted
#' @param tol mass error for deconvolution
#' @param peak_picking If TRUE, peak picking is performed.
#' @return Processed and deconvolved ESI-MS spectra.
#' @examples
#' @export

Deconvolved_ESI_High<-function(exp,maxMass,SNR.Th,max.charge=9,tol=0.015,mz1=1000,peak_picking=FALSE){

  if(peak_picking==TRUE ){


    peakDetectionCWT = getFromNamespace("peakDetectionCWT", "MassSpecWavelet")
    exp <-exp[which(exp[, 1] > mz1), ]
    Resulted<-exp
    peakInfo <- peakDetectionCWT(Resulted[,2],SNR.Th=SNR.Th,nearbyPeak=TRUE )
    majorPeakInfo = peakInfo$majorPeakInfo
    peakIndex <- majorPeakInfo$peakIndex
    betterPeakInfo <- tuneInPeakInfo(Resulted[,2], majorPeakInfo)
    mz<- Resulted[,1][betterPeakInfo$peakIndex]
    DetectedSmoothed<-cbind(mz,betterPeakInfo$peakValue)
    ordered.mz.snr<-as.matrix(DetectedSmoothed)
    ordered.mz.snr<-ordered.mz.snr[order(-ordered.mz.snr[,2]),]
    plot(DetectedSmoothed,type="h",xlim=c(3099,3101))

}else{

  ordered.mz.snr<-exp[order(-exp[,2]),]

  ordered.mz.snr<-as.matrix(ordered.mz.snr)
}
  plot(ordered.mz.snr,type="h")

  Finalresults<-data.frame()
  MONO <- c()
  Finalresultsb<-data.frame()
  ordered.transformedb<-c()
    bucle<-nrow(ordered.mz.snr)

  maxMass=maxMass
  myList<-vector("list",100)
  List_masserror<-vector("list",100)
  ListMw<-vector("list",100)
  List_deconvoluted<-vector("list",100)


  ordered.transformed_jointb<-c()
  ordered.transformed_joint<-c()

  for(KKIKIK in 1:bucle){


    if(class(ordered.mz.snr)=="matrix" ){

      if(nrow(ordered.mz.snr)>0){
        hfull<-ordered.mz.snr[1,]



    hmz<- as.numeric(hfull[1])
    hsnr<- as.numeric(hfull[2])

    if (!is.na(hmz)){

      W=1.1
      determine.charge<- subset(ordered.mz.snr,ordered.mz.snr[, 1] >  hmz-W & ordered.mz.snr[, 1] <  hmz+W)
      determine.charge
      if(nrow(  determine.charge  )>1){
        ordered.determine.charge <-determine.charge[order(-determine.charge[,1]),]
        ordered.determine.charge
        difscharges<-abs(diff(ordered.determine.charge))
        maxChargepossible<- max(1/difscharges[,1])
        maxChargepossible<- round( maxChargepossible)


        if (maxChargepossible>70){
          maxChargepossible<-max.charge

        }


        hydrogen=1.00794--0.0005485799


        b <- c()
        score_champion=0
        vectorCharge<-c(1:maxChargepossible)
        dat <- as.data.frame(matrix(ncol=5, nrow=0))
        datj <- as.data.frame(matrix(ncol=6, nrow=0))

        M3<- as.data.frame(matrix(ncol=3, nrow=0))


        Mwaa<-c()
        Mw_neutral<-c()
        Mwaa_previa<-c()
        data_champion2<-c()
        data_champion<-c()

        datleft <- as.data.frame(matrix(ncol=5, nrow=0))
        datder <- as.data.frame(matrix(ncol=5, nrow=0))

        for(i in 1:length(vectorCharge)){
          dat <- as.data.frame(matrix(ncol=5, nrow=0))
          datleft <- as.data.frame(matrix(ncol=5, nrow=0))
          datder <- as.data.frame(matrix(ncol=5, nrow=0))
          charge<-vectorCharge[i]
          m=charge*hmz
          remainingcharge<-vectorCharge[-i]
          calc.mz<-m/remainingcharge
          semi<-cbind(charge=remainingcharge,mz=calc.mz)
          full<-cbind(charge=charge,mz=hmz)
          fullx<-rbind(full,semi)
          totalscore=0
          score=0
          charge1<-fullx[1,1]*fullx[1,2]
          Mz<-as.numeric(charge1)
          data_champion3<-c()
          data_champion3_ordered <-c()

          if( Mz<=maxMass){

            W=tol
              exp.mz <- as.data.frame(matrix(ncol=0, nrow=1))
              right<-c()
              p=0
              while(nrow(exp.mz)==1 ){

                 p=p+1
                  right  <- hmz+p/charge


                   exp.mz <- subset(ordered.mz.snr,ordered.mz.snr[, 1] > right-W & ordered.mz.snr[, 1] < right+W)

                  if(nrow(exp.mz)==1){
                    score<-log(exp.mz[,2])
                    datder[p,1]<-exp.mz[,1]
                    datder[p,2]<-score
                    datder[p,3]<-exp.mz[,2]
                    datder[p,4]<-charge
                    datder[p,5]<-exp.mz[,1]* charge-charge*hydrogen
                    datder<-datder[complete.cases(datder),]
                  }

                  }

              exp.mz <- as.data.frame(matrix(ncol=0, nrow=1))
              left <-c()
              p=0
              while(nrow(exp.mz)==1 ){
                p=p+1

                  left  <- hmz-p/charge
                  exp.mz <- subset(ordered.mz.snr,ordered.mz.snr[, 1] > left-W & ordered.mz.snr[, 1] < left+W)

                  if(nrow(exp.mz)==1){
                    score<-log(exp.mz[,2])
                    datleft[p,1]<-exp.mz[,1]
                    datleft[p,2]<-score
                    datleft[p,3]<-exp.mz[,2]
                    datleft[p,4]<-charge
                    datleft[p,5]<-exp.mz[,1]* charge-charge*hydrogen
                    datleft<-datleft[complete.cases(datleft),]


                  }


              }

          dat<-c()

          dat <- rbind(datder,datleft)

            totalscore<-sum(dat[,2])


            if (totalscore>score_champion){
              score_champion=totalscore
              data_champion<-dat[,4]
              data_champion2<-dat

              Mwaa<-data_champion2[which(data_champion2[,1]==1),]
              Mw_neutral<-Mwaa[2]-hydrogen

            }

          }

          b<-c(b,totalscore)

        }



      if(!is.null(data_champion2)){
        .mz <-c()

        plot(data_champion2$V1,data_champion2$V2,type="h")

        W=1.01
      .mz <- subset(ordered.mz.snr,ordered.mz.snr[, 1] > hmz-W & ordered.mz.snr[, 1] < hmz+W)

    while(nrow(.mz)>1){
     W<-W-0.01
      .mz <- subset(ordered.mz.snr,ordered.mz.snr[, 1] > hmz-W & ordered.mz.snr[, 1] < hmz+W)

    }

        score<-mean(log(.mz[,2]))
      datj[1,1]<-.mz[,1]
      datj[1,2]<-score
      datj[1,3]<-.mz[,2]
      datj[1,4]<- data_champion[1]
      datj[1,5]<-.mz[,1]* data_champion[1]-data_champion[1]*hydrogen
      datj[1,6]<- score_champion
      data_champion3 <- rbind(datj[,-6],data_champion2)

      data_champion3
      data_champion3_ordered <-c()
      data_champion3_ordered  <- data_champion3[order(data_champion3[,1]),]

      data_champion3_ordered<- cbind(data_champion3_ordered,score_champion)


      if (!is.null(data_champion3)){

        transformed<-c()

        transformed <- exp[which(exp[, 1] >= min(data_champion3[,1]) & exp[, 1] <= max(data_champion3[,1])),]

if(class(transformed)=="data.frame"){
        transformed<-data.frame(transformed,transformed[,1]*data_champion2[1,4])

}else{

  transformed<-data.frame(transformed[1],transformed[2],as.numeric(transformed[1])*data_champion2[1,4])

colnames(transformed)<-c("mz","V2", "transformed...1....data_champion2.1..4.")
}
        ordered.transformed<-transformed[order(transformed[,3]),]

        ordered.transformedb<-rbind(ordered.transformedb,ordered.transformed)
        transformedb<-ordered.transformedb[order(ordered.transformedb[,3]),]

        plot(ordered.transformed[,3],ordered.transformed[,2],type="h")


        Finalresultsb<-rbind(Finalresultsb,data_champion3)
        Finalresults<-rbind(Finalresults,data_champion3_ordered[1,])

      }
      RR<-setdiff(ordered.mz.snr[,1],data_champion3[,1])

      }else{

        RR<-setdiff(ordered.mz.snr[,1],hmz)


      }

         }else{

        RR<-setdiff(ordered.mz.snr[,1],hmz)


      }



      ordered.transformed_jointb<- c(ordered.transformed_jointb,ordered.transformed[,3])
      ordered.transformed_joint<- c(ordered.transformed_joint,ordered.transformed[,2])
      ordered.mz.snr<-ordered.mz.snr[ordered.mz.snr[,1] %in% RR,]

    }


      }


     }


  }
  Finalresults
  Finalresultsb


if(length(Finalresults)>0){
colnames(Finalresults)<-c("Monoisotopic mass", "Score", "Intensity", "Charge","Neutral mass","Total Score")
colnames(Finalresultsb)<-c("Monoisotopic mass", "Score", "Intensity", "Charge","Neutral mass")

}

  ordered.transformed_jointFINAL<- cbind(ordered.transformed_jointb, ordered.transformed_joint)

  ordered.transformed_jointFINAL<-  ordered.transformed_jointFINAL[order(-  ordered.transformed_jointFINAL[,1]),]

plot((ordered.transformed_jointFINAL),type="h")



Totallyresults<-list(Finalresultsb,Finalresults,ordered.transformed_jointFINAL)
return(Totallyresults)
}

