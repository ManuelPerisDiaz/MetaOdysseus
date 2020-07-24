#' Pre-process and deconvolved ESI-MS spectra with simulation algorithm
#'
#' @param file Experimental MS file.
#' @param t1 Select lower m/z limit.
#' @param t2 Select higher m/z limit.
#' @param norm Normalize MS spectra. Default to TRUE.
#' @param comb Select pre-processing method.
#' @param SNR.Th Signal-to-noise ratio. Default to 3.
#' @param plotDetection Plot detected peaks after peak picking. Default to FALSE.
#' @param plotProcess Plot processed MS spectra. Default to FALSE.
#' @param plotOriginal Plot raw MS spectra. Default to TRUE.
#' @param xy MS spectra in xy file?. Default to FALSE.
#' @param dir Directory where to save plots. Default to C:.
#' @param chargemin Minimum charge state
#' @param chargemax Maximum charge state
#' @param maxMass Maximum mass deconvoluted
#' @param minMass minimum mass deconvoluted
#' @param FWHM FWHM for deconvolution
#' @param interval mass interval for deconvolution
#' @param Wa mass error for deconvolution
#' @param top Top peaks selected for deconvolution
#' @return Processed and deconvolved ESI-MS spectra.
#' @examples
#' @export

ProcessESI_Decon<-function(file,t1=300,t2=1800,norm=TRUE,comb="none",plotDetection=FALSE,plotProcess=FALSE,plotOriginal=TRUE,xy=FALSE,dir="C:"

                          , chargemin=1,chargemax=8,maxMass=8000,minMass=5000,FWHM=2,interval=5,Wa=0.5,top=100,SNR.Th=2

                           ){


  if (xy==TRUE){
    peak<-data.frame(file)

  }else{
    filex<-openMSfile(file)
    spec1<-spectra(filex,1)
    peak<-spec1
  }

  if (t1 < 0 | t2 < 0 | t1 > t2) {
    stop("No possible values for t1 and t2")
  }


  if (t2 > t1) {
    peakB <- peak[which(peak[, 1] > t1 & peak[, 1] < t2),
                  ]
  }

  smoothDWT  = getFromNamespace("smoothDWT", "MassSpecWavelet")
  sav.gol = getFromNamespace("sav.gol", "MassSpecWavelet")
  if (comb== "comb1") {
    ################ process the spectra, baseline and smooth
    ## comb1 smooth savitsky golay + baseline ALS
    v<-sav.gol(peakB[,2],fl=200)
    ycorr<-baseline.corr(v)
    Resulted<-cbind(peakB[,1],ycorr)
    # plot(Resulted,type="l
    #     ")
  }

  if (comb== "comb2") {
    ## comb2  baseline ALS +smooth savitsky golay
    ycorr<-baseline.corr(peakB[,2])
    v<-sav.gol(ycorr,fl=50)
    Resulted<-cbind(peakB[,1],v)
   # plot(Resulted,type="l
    #     ")
    }

  if (comb== "comb3") {
    ## comb3 smooth DIFSM + baseline ALS
    smoothed<-difsm(peakB[,2], lambda = 1e7)
    ycorr<-baseline.corr(smoothed)
    Resulted<-cbind(peakB[,1],ycorr)
    #  plot(Resulted,type="l
     #     ", xlim=c(1600,1610))
    }

  if (comb== "comb4") {
    ## comb4  baseline ALS +smooth DIFSM
    ycorr<-baseline.corr(peakB[,2])
    smoothed<-difsm(ycorr, lambda = 1e7)
    Resulted<-cbind(peakB[,1],smoothed)
    # plot(Resulted,type="l
    #      ")
    }
  if (comb== "comb5") {
    ## comb5  baseline ALS +smooth DWT
    ycorr<-baseline.corr(peakB[,2])
    smoothed<-smoothDWT(ycorr)
    Resulted<-cbind(peakB[,1],smoothed)
    plot(Resulted,type="l
         ")
    }
  if (comb== "none"){

    Resulted<- peakB
      plot(Resulted,type="l")
    }


  peakB<-Resulted
  colnames(peakB)<-c("V1","V2"
  )
  path_out=dir
  write.csv2(Resulted,paste(path_out,"Preprocessed_spectra.csv",sep=""))


  if (plotProcess== TRUE) {
    plot_processed(Resulted,norm=norm,dir=dir,csv=FALSE)
  }

  # tiff("profiled.tiff",width=4000, height=2500,res=500)
  # par(mar=c(5,6,4,1)+.1)
  # plot(Resulted,type="l",cex.lab=2,cex.axis=1.5,xlab = "m/z",ylab="Intensity")
  # dev.off()

  #

  ################
  #### PEAK PICKING

  peakInfo <- peakDetectionCWT(Resulted[,2],SNR.Th=SNR.Th,nearbyPeak=TRUE)
  majorPeakInfo = peakInfo$majorPeakInfo
  peakIndex <- majorPeakInfo$peakIndex
  betterPeakInfo <- tuneInPeakInfo(Resulted[,2], majorPeakInfo)


  mz<- Resulted[,1][betterPeakInfo$peakIndex]

  DetectedSmoothed<-cbind(mz,betterPeakInfo$peakValue)



  mz.snr<-data.frame(mz,betterPeakInfo$peakSNR)
  #plotPeak(Resulted[,2],betterPeakInfo$peakIndex)

  # readline(prompt="Peak detection: Press [enter] to accept it and continue")


  ##filter the top 20

  mz.snr <- mz.snr[with(mz.snr,order(-mz.snr$betterPeakInfo.peakSNR)),]

  nrow <- nrow(mz.snr)
  topB=top/100

  mz.snr <-mz.snr[1:round(topB*nrow),]

  # tiff("centroided.tiff",width=4000, height=2500,res=500)
  # par(mar=c(5,6,4,1)+.1)
  plot(mz.snr$mz,mz.snr$betterPeakInfo.peakSNR,type="h",cex.lab=2,cex.axis=1.5,xlab = "m/z",ylab="SNR")
  #dev.off()

  # readline(prompt="Peak filtering: Press [enter] to accept it and continue")

  if (plotDetection== TRUE) {
    plot_detection(DetectedSmoothed,orig=FALSE,QC2 =FALSE,QC1=FALSE,norm=norm,dir=dir)

  }

  if (plotOriginal ==TRUE){

    plot_Peak(peakB,norm=norm,dir=dir)
  }

  deconvoluted<-c()



  ordered.mz.snr<-mz.snr[order(-mz.snr$betterPeakInfo.peakSNR),]

  ordered.transformedb<-c()

  bucle<-nrow(ordered.mz.snr)
  chargemax=chargemax
  maxMass=maxMass
  minMass=minMass
  myList<-vector("list",100)
  List_masserror<-vector("list",100)
  ListMw<-vector("list",100)
  List_deconvoluted<-vector("list",100)
  List_deconvoluted_full<-vector("list",100)
  deconvoluted_Charge_LIM_lista<-vector("list",100)
  y=0
  FINALSUMA=0
  interval<-interval
  deconXX<-c()
  deconvoluted_fin_fin<-c()

  for(KKIKIK in 1:bucle){
    W=Wa
    hfull<-ordered.mz.snr[1,]
    hmz<-hfull[1]$mz
    hsnr<-hfull$betterPeakInfo.peakSNR

    if (!is.na(hmz)){

      SNRthreshold<-SNR_loop(chargemax=chargemax,maxMass=maxMass,W=W, SNRuser=1,ordered.mz.snr,hmz)

      hydrogen=1.00794--0.0005485799
      b <- c()
      score_champion=0
      score_champion2=0
      champion_consecutivecharges=0
      vectorCharge<-c(1:chargemax)
      dat <- as.data.frame(matrix(ncol=5, nrow=0))
      M3<- as.data.frame(matrix(ncol=3, nrow=0))
      data_champion2<-c()
      data_champion<-c()
      Mwaa<-c()
      Mw_neutral<-c()
      Mwaa_previa<-c()
      data_champion2_previa<-c()
      data_champion_previa<-c()

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
        dat <- as.data.frame(matrix(ncol=5, nrow=0))
        charge1<-fullx[which(fullx[, 1] ==1 ),]
        Mz<-as.numeric(charge1[2])

        if( Mz<=maxMass & Mz> minMass){

          for(i in 1:length(fullx[,1])){
            calc.mz<-as.numeric(fullx[i,2])
            CHarge<-as.numeric(fullx[i,1])

            exp.mz <- ordered.mz.snr[which(ordered.mz.snr[, 1] > calc.mz-W & ordered.mz.snr[, 1] < calc.mz+W),
                                     ]
            while(nrow(exp.mz)>1 ){
              W<-W-0.01
              exp.mz <- ordered.mz.snr[which(ordered.mz.snr[, 1] > calc.mz-W & ordered.mz.snr[, 1] < calc.mz+W),
                                       ]
            }

            if(nrow(exp.mz)==1 ){

              if( exp.mz[,2]> SNRthreshold){
                score<-log(exp.mz[,2])
                dat[i,1]<-exp.mz[,1]
                dat[i,2]<-score
                dat[i,3]<-exp.mz[,2]
                dat[i,4]<-CHarge
                dat[i,5]<-exp.mz[,1]* CHarge-CHarge*hydrogen
                dat<-dat[!duplicated(dat),]
                dat<-dat[complete.cases(dat),]
                totalscore<-totalscore+score
              }
            }

            if (totalscore>score_champion2){
              score_champion2=totalscore
              data_champion_previa<-fullx
              data_champion2_previa<-dat
              Mwaa_previa<-data_champion_previa[which(data_champion_previa[,1]==1),]
              Mw_neutral_previa<-Mwaa_previa[2]-hydrogen
            }

#here we filter for maximum difference 2
            if (length(data_champion2_previa)>0) {
                if (nrow(data_champion2_previa)>1){

                   numbers=data_champion2_previa$V4[order(data_champion2_previa$V4)]

                   jumps2<-abs(diff(numbers))

                   d1<-which(jumps2 < 2)
                   indexes<-numbers[c(d1, d1[length(d1)]+1)]
                   data_champion2_previa<-  subset(data_champion2_previa,data_champion2_previa$V4 %in%indexes )

                    consecutivecharges<-sum((jumps2==1))

                    if(consecutivecharges>champion_consecutivecharges){
                      champion_consecutivecharges=consecutivecharges
                      score_champion=score_champion2
                      data_champion<-data_champion_previa
                      data_champion2<-data_champion2_previa
                      Mwaa<-data_champion_previa[which(data_champion_previa[,1]==1),]
                      Mw_neutral<-Mwaa[2]-hydrogen

                    }


                }
}

          }

          b<-c(b,score_champion)

        }

      }


      ### we have the first peak series, for the most intensitve peak in the spectrum, now
      ## we are simulating the spectrum
      #so transform signals in the peak to the zero-charge mass domain, select sample every X mass
      ordered.transformed<-c()
      transformed<-c()
      Finalresults<-c()
      dat<-data_champion2
      if(length(dat)>0){
      if (nrow(dat)>0  ){
        m=mean(dat[,5])
        ## we need to exctract for every peak, for every row of dat, the part of the mass spectrum
        ## the broadness is determined by FWHM, leave at 20.
        #interval<-interval
        transformed<-c()


        transformedC <- peakB[which(peakB[, 1] > hmz-interval & peakB[, 1] < hmz+interval),]
        transformed<-data.frame(transformedC,transformedC[,1]*data_champion2[1,4]-data_champion2[1,4]*hydrogen)


        ordered.transformed<-transformed[order(transformed[,3]),]

        ordered.transformedb<-rbind(ordered.transformedb,ordered.transformed)
        transformedb<-ordered.transformedb[order(ordered.transformedb[,3]),]
        transformedbdeconc<-ordered.transformedb[order(ordered.transformedb[,3]),]

        #  plot(ordered.transformed[,3],ordered.transformed[,2],type="l",main="13")
        # plot(transformedb[,3],transformedb[,2],type="l",main="Deconvoluted M")
        # plot(transformedb[,3],transformedb[,2],type="h",main="Deconvolved M")

          ##  tiff("deconvoluted.tiff",width=4000, height=2500,res=500)
          ##   par(mar=c(5,6,4,1)+.1)
          ##   plot(transformedb[,3],transformedb[,2],type="l",cex.lab=2,cex.axis=1.5,xlab = "Mass (Da)",ylab="Intensity")
          ##   dev.off()

        ## ## ### denote these peaks as processed and eliminate them from ordered


        dat<-dat[!duplicated(dat),]


        Finalresults<-rbind(Finalresults,dat)
      }


      RR<-setdiff(ordered.mz.snr[,1],dat[,1])

      ordered.mz.snr<-ordered.mz.snr[ordered.mz.snr[,1] %in% RR,]

      ####Peak list, now deconvolute isotopes

      FInaldupl<-Finalresults[!duplicated(Finalresults$V1),]

      Finbbb<-data.frame()

      for (L in 1:nrow(FInaldupl)){
        Finb<-c()
        Finbb<-data.frame()
        FIn<-c()
        M<-FInaldupl[L,1]
        ChargeM<-FInaldupl[L,4]
        MassFInM<-FInaldupl[L,5]

        ## now we  determine the average mass from the isotopic distribution
        for (W in 1:200){
          FIn<- FInaldupl[which(FInaldupl[, 1] > (M-W/ChargeM) & FInaldupl[, 1] < (M+W/ChargeM)),]
          if (nrow(FIn)>0){
            for (j in 1:nrow(FIn)){
              if (FIn[j,4]==ChargeM ){

                if( (abs(FIn[j,5]-MassFInM))< 5 ){
                  Finb<-rbind(Finb,FIn[j,])
                  Finb<-Finb[!duplicated(Finb[,1]),]
                  Finbb<-data.frame(mean(Finb[,1]),sum(Finb[,2]),sum(Finb[,3]),mean(Finb[,4]),mean(Finb[,5]))

                } }}} }

        Finbbb<-rbind(Finbbb,Finbb)


      }

      ## we need to exctract for every peak, for every row of dat, the part of the mass spectrum
      ## the broadness is determined by FWHM, leave at 20.

      Data.FRAMEE<-c()
      mu<-c()
      sd<-c()
      peakheight<-c()
      charge<-c()
      for (i in 1:length(Finbbb[,1])){

        transformed <- peakB[which(peakB[, 1] > Finbbb[i,1]-FWHM & peakB[, 1] < Finbbb[i,1]+FWHM),]

        while (length(transformed[,2]) <5){
          FWHM<-FWHM+0.5
          transformed <- peakB[which(peakB[, 1] > Finbbb[i,1]-FWHM & peakB[, 1] < Finbbb[i,1]+FWHM),]

           print("Increase FWHM, increased +0.5")

        }

        # plot(transformed,type="l",main="Extracted charge from m/z domain with FWHM window")

        GaussianFitted<-NULL
        try(silent=TRUE,GaussianFitted<-FitPeakToGaussian(transformed[,1],transformed[,2],Finbbb[i,1],doPlot = TRUE,      xTitle='m/z', yTitle='Intensity')
        )

        if(!is.null(GaussianFitted)){
          mu[i]<-GaussianFitted$coefficients[1]
          sd[i]<- GaussianFitted$coefficients[2]
          peakheight[i]<- GaussianFitted$coefficients[3]
          charge[i]<-Finbbb[i,4]
        }else{

          pk.pos <-findpeaks(transformed[,2])

          if (length(pk.pos)>0){
          pks<-fitpeaks(transformed[,2],pk.pos)

          if(any(!is.na(pks))){
          mas<-which.max(pks[,4])

          mu[i]<-transformed[pk.pos[mas],1]
          sd[i]<- pks[2]
          peakheight[i]<- pks[4]
          charge[i]<-Finbbb[i,4]
          }

          }
          #  plot(transformed,type="l",main="10 fitted second way")
          # apply(pks, 1,
          #       function(pkmodel) {
          #         lines(transformed[,1],
          #               dnorm(1:length(transformed[,1]), pkmodel["rt"], pkmodel["sd"]) *
          #                pkmodel["area"],
          #              col = 2)
          #         invisible()
          # })

        }
      }


        Data.FRAMEE<-data.frame(mu=mu,sd=sd,peakheight=peakheight,charge)

      if(length(Data.FRAMEE)>0){
      envelope<-data.frame(charge,peakheight)
      envelope<-envelope[order(-envelope[,1]),]
      x<-envelope[,1]
      y<-envelope[,2]

      #  tiff("charge.tiff",width=4000, height=2500,res=500)
      # par(mar=c(5,6,4,1)+.1)
       #plot(x,y,type="h"
        #  ,cex.lab=2,cex.axis=1.5,xlab = "Charge state",ylab="Intensity")
       # dev.off()

      Mwaverage<-mean((Data.FRAMEE[,1]*Data.FRAMEE[,4])-(Data.FRAMEE[,4]*hydrogen))

      masserror<-sd((Data.FRAMEE[,1]*Data.FRAMEE[,4])-(Data.FRAMEE[,4]*hydrogen))

      mu <-(Mwaverage + (Data.FRAMEE[,4]*hydrogen))/Data.FRAMEE[,4]

      x<-Resulted[,1]
      yc2<-as.data.frame(matrix(ncol=6, nrow=length(x)))
      deconvoluted_Charge<-as.data.frame(matrix(ncol=6, nrow=length(x)))
      deconvoluted_M<-as.data.frame(matrix(ncol=6, nrow=length(x)))
      yc2[,1]<-x

      for(i in 1:length(mu)){
        yc2[,i+1] <- peakheight[i]*exp(-0.5*(x-mu[i])^2/sd[i]^2)

       # tiff("charge.tiff",width=4000, height=2500,res=500)
        # par(mar=c(5,6,4,1)+.1)
         plot(yc2[,1],yc2[,i+1],type="l",
              cex.lab=2,cex.axis=1.5,xlab = "m/z",ylab="Intensity")
         #dev.off()

        deconvoluted_Charge[,i]<-Data.FRAMEE[i,4]
        deconvoluted_M[,i]<-x*Data.FRAMEE[i,4]
      }

      deconX<-c()
      deconvoluted_fin<-c()
      deconvoluted<-c()
       for(i in 1:length(mu)){
        if( is.finite(yc2[,i+1]) ){
            plot(deconvoluted_M[,i],yc2[,i+1],type="l",main="Deconvoluted Simulated individual m/z")
       deconvoluted<-c(yc2[,i+1],deconvoluted_M[,i])

        }
         deconX<-c(deconX,deconvoluted_M[,i])
         deconvoluted_fin<-c(deconvoluted_fin,yc2[,i+1])
         }

      Zdec<-data.frame(deconX,deconvoluted_fin)


      #tiff("deconv.tiff",width=4000, height=2500,res=500)
      #par(mar=c(5,6,4,1)+.1)

       #plot(deconX,deconvoluted_fin,type="l",xlim=c(6000,6100),cex.lab=2,cex.axis=1.5,xlab = "Mass(Da)",ylab="Intensity")
       #dev.off()
       #library(plyr)
        poiu<-ddply(Zdec, "deconX", numcolwise(sum))



      fil<-yc2[colSums(!is.na(yc2))>0]
      if (ncol(fil)>2){
        suma<-rowSums(fil[,2:length(fil[,])])
        Fullserie<-cbind(fil,suma)
      }
      if (ncol(fil)==2){
        suma<-fil
        Fullserie<- suma
        names(Fullserie)[2] <- "suma"

      }


      deconvoluted_Charge_LIM<-deconvoluted_Charge[colSums(!is.na(deconvoluted_Charge))>0]


      ## extract rows where column 2 is > que 0 for the Fullserie[,5], these are the m/z should be
      ## substracted from the peak[,1] vector

   # tiff("sim_exp_one.tiff",width=4000, height=2500,res=500)
      #   par(mar=c(5,6,4,1)+.1)
     # plot(col="red",Fullserie[,1],Fullserie$suma,type="l",cex.lab=2,cex.axis=1.5,xlab = "m/z",ylab="Intensity")


      # plot(col="red",peak,type="l",cex.lab=2,cex.axis=1.5,xlim=c(960,1800),xlab = "m/z",ylab="Intensity")
      #lines(Fullserie[,1],Fullserie$suma,col="black")

      #dev.off()

      Fullserie_extracted<- Fullserie[which(Fullserie$suma > 0.000001),]
     # plot(Fullserie_extracted[,1],Fullserie_extracted$suma,type="h",main="4")
      RRdif<-setdiff(peakB[,1],Fullserie_extracted[,1])
      peakBH<-peakB[peakB[,1] %in% RR,]
      #plot(peakBH,type="h")

      myList[[KKIKIK]]<-Fullserie
      ListMw[[KKIKIK]]<-Mwaverage
      List_masserror[[KKIKIK]]<-masserror
      List_deconvoluted[[KKIKIK]]<-Finbbb

      List_deconvoluted_full[[KKIKIK]]<-deconvoluted_M

      deconvoluted_Charge_LIM_lista[[KKIKIK]]<-yc2

      #substrcc<-peakB[,2]-Fullserie$suma

      Combined<-data.frame(peakB,Fullserie$suma)

      resta<-Combined[,2]-Combined$Fullserie.suma

      peakB<-data.frame(cbind(Combined$V1,resta))

       peakB[peakB<0] <- 0

     names(peakB)[1]<-"V1"
      names(peakB)[2]<-"V2"
     # plot(peakB,type="l",main="1")

      }

      }else{
        ##remove this element from the list
        metocas<-setdiff(ordered.mz.snr[,1],hfull[,1])
        ordered.mz.snr<-ordered.mz.snr[ordered.mz.snr[,1] %in% metocas,]
      }

    }

    deconXX<-c(deconXX,deconX)
    deconvoluted_fin_fin<-c(deconvoluted_fin_fin,deconvoluted_fin)
    finjl<-data.frame(deconXX, deconvoluted_fin_fin)
    hfull<-ordered.mz.snr[1,]
    hmz<-hfull[1]$mz

    if (is.na(hmz)){
      return(list(myList,ListMw,List_masserror,List_deconvoluted,List_deconvoluted_full,deconvoluted_Charge_LIM_lista,Resulted,transformedbdeconc,deconvoluted, finjl))
    }

  }


}

