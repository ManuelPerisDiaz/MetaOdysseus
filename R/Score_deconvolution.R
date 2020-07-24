#' Scoring function for peptide mass fingerprint MALDI-MS experiment.
#'
#' @param ESIprocessed Processed ESI-MS spectra.
#' @param intervaLo mass interval for deconvolution
#' @param intervaLo2 mass interval for deconvolution
#' @param loess If TRUE, loess algorithm is used.
#' @param span span controlling the loess
#' @param mz1 Left mass
#' @param mz2 Right mass
#' @examples
#' @return Annotation score for PMF MALDI-MS.
#' @export

Score_deconvolution<-function(ESIprocessed, intervaLo=3000,
                              intervaLo2=2000,loess=FALSE,span=0.002, mz1=NULL,mz2=NULL){

  if (mz1 < 0 | mz2 < 0 | mz1 > mz2) {
    stop("No possible values for t1 and t2")
  }

  options(warn=-1)
  lista<-ESIprocessed
  hydrogen=1.00794--0.0005485799

  Uscore_lista<-vector("list",100)
  DScoretodos_lista<-vector("list",100)
  totalFscore_lista<-vector("list",100)
  CCScoretodos_lista<-vector("list",100)
  MScore_weighted_lista<-vector("list",100)

  Mw_lista<-vector("list",100)


### ZERO-CHARGE MASS DISTRIBUTION SPECTRUM
## Obtain for each mass, the deconvoluted total peak and the FWHM

if(loess=="TRUE"){

  x<-lista[[10]][,1]
  y<-lista[[10]][,2]
  dat<-data.frame(x,y)

  poiu<-ddply(dat, "x", numcolwise(mean))
  poiu<-aggregate(dat$y, by=list(dat$x), FUN=mean)

  plot(poiu,type="h",xlim=c(7000,7200))

  loessData <- data.frame(
    x = poiu$Group.1,
    y = predict(loess(poiu$x~poiu$Group.1, dat, span = span)),
    method = "loess()"
  )


  M<-loessData$x
  Y<-loessData$y
  Deconvoluted<-data.frame(M,Y)
  DEb<-data.frame(loessData$x,loessData$y)
  plot_Peak(DEb,ID="Deconvoluted_Loes", norm = TRUE,t1 = mz1, t2 = mz2,axis=c("Mass(Da)"))

  plot(loessData$x,loessData$y,type="l",xlim=c(7000,7500))



  }else{
    x<-lista[[8]][,3]
    y<-lista[[8]][,2]


    dat<-data.frame(x,y)

    approxData <- data.frame(
      with(dat,
           approx(x, y, xout = seq(min(x), max(x), by = 0.0001), method = "linear")
      ),
      method = "approx()"
    )

M<-approxData$x
Y<-approxData$y
Deconvoluted<-data.frame(M,Y)

DE<-data.frame(approxData$x,approxData$y)
plot_Peak(DE,ID="Deconvoluted",norm = TRUE, t1 = mz1, t2 = mz2,axis=c("Mass(Da)"))
plot(approxData$x,approxData$y,type="l",xlim=c(7100,7500))

}
Deconvolvedmass<-c()
for(i in 1:length(lista[[2]])){
  Deconvolvedmass<-c(Deconvolvedmass,lista[[2]][[i]])

}


indexesDeconvolvedmass<-match(Deconvolvedmass,lista[[2]] )
mass<-lista[[5]]
intensi<-lista[[6]]
a<-c()
FWHM <- c()
intervaLo<-intervaLo
intervaLo2<-intervaLo2
UScore_weighted<-c()
SUMAAAAtotal<-vector("list",100)
totalFscore<-c()
indexesDeconvolvedmass_b<-c()
indexesDeconvolvedmass<-indexesDeconvolvedmass[which(!is.na(indexesDeconvolvedmass))]

for (j in indexesDeconvolvedmass){

  Fscore<-c()
  z<-c()
  m<-c()

  Mw<-mean(lista[[4]][[j]][,5])
  Mw_lista[[j]] <-Mw
  M<-mass[[j]]
  In<-intensi[[j]]
  M <- M[,colSums(is.na(M))<nrow(M)]
  In <- In[,colSums(is.na(In))<nrow(In)]
  In<-In[,-1]

  Xtotal<-c()
  Xtot<-c()
  Dta <-c()
  Colyy<-c()
  data.<-c()
  if (is.numeric(M)){
    M<-data.frame(M)

  }
  if (is.numeric(In)){
    In<-data.frame(In)

  }


  for(i in 1:ncol(In)){


    y<-c(In[,i])
    z<-c(z,lista[[4]][[j]][i,4])
    m<-c(m,lista[[4]][[j]][i,1])

    x<-c(M[,i])-hydrogen*z[i]

    XX<-data.frame(x,y)


    Coly<-subset(XX$y,XX$x>Mw-intervaLo & XX$x< Mw+intervaLo)
    Colx<-subset(XX$x,XX$x>Mw-intervaLo & XX$x< Mw+intervaLo)
    Colyy<-cbind(Colyy,Coly)


    plot(x,y,type="l",main="Zero-charge mass distribution Individual charge")

  }

  Dta <-cbind(Colx,Colyy )
  SUMAAAA<-c()
   SUMAAAA <- rowSums(Colyy)
   SUMAAAAtotal[[j]]<-SUMAAAA
  Deconvolutedb<-Deconvoluted[which(Deconvoluted[, 1] > Mw -intervaLo2 &  Deconvoluted[, 1] <  Mw +intervaLo2),]

  while (nrow(Deconvolutedb) <10){

    intervaLo2<-intervaLo2+5
    Deconvolutedb<-Deconvoluted[which(Deconvoluted[, 1] > Mw -intervaLo2 &  Deconvoluted[, 1] <  Mw +intervaLo2),]
    print("Increase intervaLo2, increased +5")

  }

   plot(Deconvolutedb[,1],Deconvolutedb[,2],type="l",main="Zero-charge mass distribution for particular M")

   y<-Deconvolutedb[,2]
   x<-Deconvolutedb[,1]

   xmax <- x[y==max(y)]


  x1 <- x[x < xmax][which.min(abs(y[x < xmax]-max(y)/2))]
  x2 <- x[x > xmax][which.min(abs(y[x > xmax]-max(y)/2))]


  if(length(x1)>0 & length(x2)>0){
  points(c(x1, x2), c(y[x==x1][1], y[x==x2][1]), col="red")}

  if (length(x1)==0){

    Mw<-mean(lista[[4]][[j]][,5])

    Deconvolutedb<-Deconvoluted[which(Deconvoluted[, 1] > Mw -intervaLo2 &  Deconvoluted[, 1] <  Mw +intervaLo2),]



    plot(Deconvolutedb[,1],Deconvolutedb[,2],type="l",cex.lab=2,cex.axis=1.5,xlab = "Mass (Da)",ylab="Intensity",main="Zero-charge mass distribution for particular M")


    y<-Deconvolutedb[,2]
    x<-Deconvolutedb[,1]

    xmax <- x[y==max(y)]

    x1 <- x[x < xmax][which.min(abs(y[x < xmax]-max(y)/2))]
    x2 <- x[x > xmax][which.min(abs(y[x > xmax]-max(y)/2))]

    if(length(x1)>0 & length(x2)>0){
    points(c(x1, x2), c(y[x==x1][1], y[x==x2][1]), col="red",cex=5,pch=3)}

  }


  FWHM <- x2-x1
  A<-(x2-Mw)/FWHM

  B<-(Mw-x1)/FWHM


  if (length(FWHM)>0){


  ##fscore
  MaxFWHM<-max(A, B)
  if( MaxFWHM >0.75   ){
    if (A==MaxFWHM){

      Ymin<-y[x==x2]
      Fscore<-abs((Ymin-0.5*max(y))/(0.5*max(y)))

    }else{

      Ymin<-y[x==x1]

      Fscore<-abs((Ymin-0.5*max(y))/(0.5*max(y)))
    }
  }


  ## Detect peaks by local maxima with second derivate
  # isolate x values between the fwhm range of the max


  data.<-data.frame(x,y)
  x222<- subset(data.,data.[,1]>xmax-FWHM &data.[,1]< xmax+FWHM)

  Local_maxima<-which(diff(sign(diff(x222$y)))==-2)+1

  if (length( Local_maxima)>1) {

     deltaInt<- min(abs(diff(x222[Local_maxima,]$y)))


    if (deltaInt>0.5*abs(max(x222[Local_maxima,]$y))){


      Ymin<- deltaInt
      if(length(Fscore)==0){
        Fscore=1
      }

      Intensity_peaks <-abs(max(x222[Local_maxima,]$y))
      Fscore<-Fscore*(Ymin-0.5*Intensity_peaks)/(0.5*Intensity_peaks)



      }



  }

  if(length(Fscore)==0){
    Fscore=1
  }

  totalFscore_lista[[j]] <-Fscore
  totalFscore<-c(totalFscore,Fscore)

  a<-c(a,FWHM)
  indexesDeconvolvedmass_b<-c(indexesDeconvolvedmass_b,j)

  }else{

    indexesDeconvolvedmass <-indexesDeconvolvedmass[-j]

  }

}


##############################   UScore_weighted

##calculate fwhm from deconvoluted non zero peaks


mass<-lista[[5]]
intensi<-lista[[6]]
counter=1
UScore_weighted<-c()

for (j in indexesDeconvolvedmass_b){


  FWHM<-a[counter]


  z<-c()
  m<-c()
  Mw<- lista[[2]][[j]]
  M<-mass[[j]]
  In<-intensi[[j]]
  M <- M[,colSums(is.na(M))<nrow(M)]
  In <- In[,colSums(is.na(In))<nrow(In)]

  UScoretodos<-c()
  Ysquaretodos<-c()
  Ysquare<-c()

  In<-In[,-1]
  if (is.numeric(M)){
    M<-data.frame(M)

  }
  if (is.numeric(In)){
    In<-data.frame(In)

  }


  for (i in 1:ncol(In)){
    Y<-c()
    X<-c()
    UScore<-c()


    y<-c(In[,i])
    z<-c(z,lista[[4]][[j]][i,4])
    m<-c(m,lista[[4]][[j]][i,1])
    x<-c(M[,i])-hydrogen*z[i]

    scaledFWHM<-FWHM/z[i]

    individualdeconvo<-c()
    individualdeconvo<-cbind(lista[[1]][[j]][1],lista[[1]][[j]][[i+1]])

    plot(individualdeconvo[,1],individualdeconvo[,2] ,type="l")

    Y<-subset(individualdeconvo,individualdeconvo[,1]>m[i]-scaledFWHM & individualdeconvo[,1]<m[i]+scaledFWHM)
    X<-subset(lista[[7]],lista[[7]][,1]>m[i]-scaledFWHM & lista[[7]][,1]<m[i]+scaledFWHM)


    plot(X,cex.lab=2,cex.axis=1.5,xlab = "m/z",ylab="Intensity")
    lines(Y,col="red",lwd=5)

    UScore<-1-rae(X[,2],Y[,2])
    Ysquare<-c(Ysquare,sum(Y[,2])^2)
    Ysquaretodos<-c(Ysquaretodos,Ysquare)
    UScoretodos<-c(UScoretodos,UScore)

  }

  UScore_weighted<-c(UScore_weighted,weighted.mean(UScoretodos,Ysquare))
  Uscore_lista[[j]]<- weighted.mean(UScoretodos,Ysquare)

  counter=counter+1
}


##############################   MScore_weighted and CSScore
counter=1

MScore_weighted<-c()
CCScoretodos<-c()
for (j in indexesDeconvolvedmass_b){
  z<-c()
  m<-c()
  FWHM<-a[counter]

  Mw<-mean(lista[[4]][[j]][,5])

  M<-mass[[j]]
  In<-intensi[[j]]
  M <- M[,colSums(is.na(M))<nrow(M)]
  In <- In[,colSums(is.na(In))<nrow(In)]
  In<-In[,-1]
  if (is.numeric(M)){
    M<-data.frame(M)

  }
  if (is.numeric(In)){
    In<-data.frame(In)

  }
  Xtotal<-c()
  Xtot<-c()
  Dta <-c()
  Colyy<-c()
  Coly<-c()
  for(i in 1:ncol(In)){

    y<-c(In[,i])
    z<-c(z,lista[[4]][[j]][i,4])
    x<-c(M[,i])-hydrogen*z[i]
    m<-c(m,lista[[4]][[j]][i,1])
    XX<-data.frame(x,y)

    Coly<-subset(XX$y,XX$x>Mw-  FWHM & XX$x< Mw+  FWHM)
    Colx<-subset(XX$x,XX$x>Mw-  FWHM & XX$x< Mw+  FWHM)
    Colyy<-cbind(Colyy,Coly)



    plot(x,y,type="l",main="Zero-charge mass distribution Individual charge (z=4)",xlim=c(6020,6060),cex.lab=2,cex.axis=1.5,xlab = "Mass (Da)",ylab="Intensity")

  }

  Dta <-cbind(Colx,Colyy )

  Deconvolutedb<-Deconvoluted[which(Deconvoluted[, 1] > Mw -  FWHM &  Deconvoluted[, 1] <  Mw +  FWHM),]

  if (nrow(Deconvolutedb)>0){
   plot(Deconvolutedb[,1],Deconvolutedb[,2],type="l",main="Zero-charge mass distribution for particular M")
  X<- Deconvolutedb[,2]
  x<-Deconvolutedb[,1]

  UScoretodos<-c()
  Ysquaretodos<-c()
  Ysquare<-c()
  Zsum<-c()
  for (i in 1:ncol(Colyy)){
    Y<-Colyy[,i]
    output_x_vals <- seq(min(Deconvolutedb[,1]),max(Deconvolutedb[,1]),length.out=length((Deconvolutedb[,1])))
    approxDataY <- data.frame(
      with(dat,
           approx(Dta[,1], Y, xout=output_x_vals, rule=1,method = "linear")
      ),
      method = "approx()"
    )


    plot(approxDataY$x,  approxDataY$y,type="l")

    Approl<-data.frame(approxDataY$x,  approxDataY$y,X)
    Approl<-Approl[complete.cases( Approl), ]


    normX<-normalize( Approl[,3])
    normY<-normalize(Approl[,2])

    plot(Approl[,1],normX,type="l",cex.lab=2,cex.axis=1.5,xlab = "Mass (Da)",ylab="Intensity")

    lines(Approl[,1],normY,col="red",lwd=5)

    MScore_n<-1-rae( normX,normY)


    Ysquare<-c(Ysquare,sum(Y)^2)
    Ysquaretodos<-c(Ysquaretodos,Ysquare)

    UScoretodos<-c(UScoretodos,MScore_n)


    ## for ccs score
    Z<-sum(Colyy[,i])
    Zsum<-c(Zsum,Z)


  }

  }else{

    UScoretodos<-c()
    Ysquaretodos<-c()
    Ysquare<-c()
    Zsum<-c()

    for (i in 1:ncol(Colyy)){
      Y<-Colyy[,i]
      Ysquare<-c(Ysquare,sum(Y)^2)
      Ysquaretodos<-c(Ysquaretodos,Ysquare)
      Z<-sum(Colyy[,i])
      Zsum<-c(Zsum,Z)
      MScore_n<-0
      UScoretodos<-c(UScoretodos,MScore_n)
       }


     }


  MScore_weighted<-c(MScore_weighted,weighted.mean(UScoretodos,Ysquare))
  MScore_weighted_lista[[j]]<-weighted.mean(UScoretodos,Ysquare)

  counter=counter+1

  ## for ccs score
  Zsum_carga<-data.frame(z,Zsum)

  Max.Zsum_carga<-Zsum_carga[which.max(Zsum_carga[,2]),]
  Zmax<-as.numeric(Max.Zsum_carga[1])



  plot(Zsum_carga,type="h",cex.lab=2,cex.axis=1.5,xlab = "Charge",ylab="Intensity")


  left<- data.frame(subset(Zsum_carga,Zsum_carga[,1]<Zmax))
  right<-subset(Zsum_carga,Zsum_carga[,1]>Zmax)

  left<-left[order(-left[,1]),]
  right<-right[order(right[,1]),]

  B<-0
  added<-0
  min<-c()
  if(nrow(left)>0){
  for (i in 1:nrow(left)){
    left_extracted<-left[i,]
    left_extracted_2<-left[i+1,]
    left_extracted_3<-left[i+2,]


    if(left_extracted[,2]>as.numeric(Max.Zsum_carga[2])){

      added<-left_extracted[,2]-as.numeric(Max.Zsum_carga[2])
      B<-B+added
    }


    if (any(!is.na(left_extracted_3))){
      if(left_extracted_3[,2]>left_extracted[,2]){

        min<-left_extracted_3[,2]
      }else{
        min<-left_extracted[,2]

      }}
    if(length(min)==0){  min<-left_extracted[,2]}
    if (any(!is.na(left_extracted_2))){
      if(left_extracted_2[,2]>left_extracted[,2]){

        added<-left_extracted_2[,2]- min
        B<-B+added
      }}


  }}
  added<-0
  min<-c()

  if(nrow(right)>0){
  for (i in 1:nrow(right)){
    right_extracted<-right[i,]
    right_extracted_2<-right[i+1,]
    right_extracted_3<-right[i+2,]


    if (any(!is.na(right_extracted))){
      if(right_extracted[,2]>as.numeric(Max.Zsum_carga[2])){

        added<-right_extracted[,2]-as.numeric(Max.Zsum_carga[2])
        B<-B+added
      }
    }

    if (any(!is.na(right_extracted_3))){
      if(right_extracted_3[,2]>right_extracted[,2]){

        min<-right_extracted_3[,2]
      }else{
        min<-right_extracted[,2]

      }}
    if(length(min)==0){  min<-right_extracted[,2]}

    if (any(!is.na(right_extracted_2))){
      if(right_extracted_2[,2]>right_extracted[,2]){

        added<-right_extracted_2[,2]- min
        B<-B+added
      }}

  }
  }

  CCScore=1-(B/sum(Zsum))

  CCScoretodos<-c( CCScoretodos,CCScore)

  CCScoretodos_lista[[j]]<-CCScore
}


##############################  DScore
DScore<-c()
DScoretodos<-c()
for (i in 1:length(totalFscore)){


  DScore[i] <- totalFscore[i]*
    CCScoretodos[i]*
    MScore_weighted[i]*
    UScore_weighted[i]

  DScoretodos<-c(DScoretodos,DScore[i])



}
i=0
i=1
for (j in indexesDeconvolvedmass_b){
  DScoretodos_lista[[j]]<- DScoretodos[i]

  i=i+1
}



##############################  UniScore
Suma<-c()
Total<-c()
counter=1
for (j in indexesDeconvolvedmass_b){

  if(!is.null(SUMAAAAtotal[[j]])){


   Product<- DScoretodos[counter]*(SUMAAAAtotal[[j]])^2
   Total<-c(Total,Product)
    Suma<-c(Suma,(SUMAAAAtotal[[j]])^2)

    counter=counter+1
  }


}



lista<-ESIprocessed

y=0

for (i in 1:length(lista[[1]])){

  if(!is.null(lista[[1]][[i]]$suma)){
    y=y+lista[[1]][[i]]$suma
   plot(y,type="l")

  }

}

y<-baseline.corr(y)

 tiff("simulated_todos.tiff",width=4000, height=2500,res=500)
 par(mar=c(5,6,4,1)+.1)


plot(lista[[7]][,1],lista[[7]][,2],type="l",cex.lab=2,cex.axis=1.5,xlab = "m/z",ylab="Intensity")
lines(lista[[7]][,1],y,col="#f8766d")

dev.off()

lines(lista[[7]][,1],lista[[1]][[2]]$suma,col="#f8766d")


rsq <- function (x, y) cor(x, y) ^ 2

R2<-rsq(y, lista[[7]][,2])


UniScore= R2*( sum(Total)/sum(Suma))

my_list<-list(Mw_lista,UniScore,R2,Uscore_lista,DScoretodos_lista,totalFscore_lista,
               CCScoretodos_lista,MScore_weighted_lista)

names(my_list)<-c("Molecular_weight","Uniscore","R2_coeff","Uscore","DScore","Fscore"
                  ,"CCScore","MScore"
                  )
return( my_list  )
}
