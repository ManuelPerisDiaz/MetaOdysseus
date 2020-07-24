#' Annotate ESI-MS/MS spectra
#'
#' @param Deconv_CysMpro Deconvolved ESI-MS spectra.
#' @param theory Theoretical ESI-MS/MS library
#' @param mzw Mass-to-charge window (Da) for the annotation. Default to 1.
#' @param b Intensity threshold
#' @param plot If TRUE, plot the results
#' @param supress Control the plot options
#' @param xlim Select m/z range
#' @param label Label the results
#' @examples
#' @return Peaks annotated in the MALDI-MS spectra.
#' @export
#'
targetMSMS<-function(Deconv_CysMpro,theory,mzw = 1, b = 5,plot=FALSE,
                      supress = FALSE,xlim = c(100,2500),label = "")
{

 theory<- data.frame(theory)



 expt<-data.frame(mz=Deconv_CysMpro[[2]][1],intensity=Deconv_CysMpro[[2]][3])
 expt$normalized <- ((expt[, 2]-min(expt[, 2]))/(max(expt[, 2])-min(expt[, 2])))*100
 expt<-cbind(expt,carga=Deconv_CysMpro[[2]][4],Totalscore=Deconv_CysMpro[[2]][6])
 colnames(expt)<-c("mz","intensity","normalized","Charge","Total Score")


  expt <- subset(expt, expt$normalized >= b)

  m <- matrix(0, ncol =
                14)
  identifications_fin<-data.frame(m)
  colnames(identifications_fin)<- c("expt_mz","expt_intensity","ms1seq","ms1z1", "ms1z2","ms1z3","ms1z4","ms1z5","ms1z6", "ms2seq", "ms2type", "ms2mz","Nmetal","score")



  matches <- vector("list")
  for (i in 1:nrow(expt)) {
    tmp_matches <- data.frame(NULL)
    tmp_matchesb <- data.frame(NULL)
    tmp_matchesb <- theory[expt$mz[i] >= theory$ms2mz - mzw &
                            expt$mz[i] <= theory$ms2mz + mzw, ]

    if(nrow(tmp_matchesb)>0)
    for(u in 1:nrow(tmp_matchesb)){
    splited<-strsplit(as.character(tmp_matchesb[u,]$ms2type)[1],
                      split = "")[[1]]

    annotatedcharge<-splited[6]
    tmp_matches <- tmp_matchesb[u,][expt$Charge[i] == annotatedcharge, ]

    num_tmp_matches <- nrow(tmp_matches)
    expt_mz <- rep(expt$mz[i], times = num_tmp_matches)
    expt_intensity <- rep(expt$normalized[i], times = num_tmp_matches)
    score <- rep(expt$`Total Score`[i], times = num_tmp_matches)
    matches[[u]] <- data.frame(expt_mz, expt_intensity ,score ,tmp_matches)


    }

    identifications <- as.data.frame(do.call("rbind", matches))
    identifications_fin<- rbind(identifications_fin,identifications)

  }

  identifications <-  identifications_fin

   if (nrow(identifications) == 0){
     aaa<-c()}else{

  aa<-identifications[order(identifications$ms2seq,-abs(identifications$expt_intensity)),]

  aaa<-identifications[!duplicated(aa$ms2seq),]

  if (nrow(aaa) == 0){
    aaa<-c()

  }else{

  row.names(aaa) <- 1:nrow(aaa)
  aaa$error <- round(aaa$expt_mz -
                       aaa$ms2mz, digits = 3)
  num_identifications <- nrow(aaa)



  if (plot==TRUE){

  getLocation <- function(type) {
    tmp <- strsplit(as.character(type), split = "[[:punct:]]")[[1]][2]
    tmp <- strsplit(tmp, split = "[[:alpha:]]")[[1]][2]
    return(as.numeric(tmp))
  }
  location <- sapply(aaa$ms2type, getLocation)
  getSeries <- function(type) {
    tmp <- strsplit(as.character(type), split = "[[:punct:]]")[[1]][2]
    tmp <- strsplit(tmp, split = "[[:digit:]]")[[1]][1]
    return(tmp)
  }
  series <- sapply(aaa$ms2type, getSeries)
  color <- sapply(series, function(x) ifelse(x == "b" | x ==
                                               "c", "black", "black"))

  par(mar=c(4,7,2,1))
   plot(expt$mz, expt$normalized, type = "l", bg='black',cex.lab=2.5,col="red",cex.axis=2, ylim = c(0,150), xlab = "m/z", ylab = "", yaxs = "i",
       yaxt = "n")
   myTicks<-c(0.00, 0.50, 1.00)
   axis(2,at = c(0.00, 50.00, 100.00), labels = sprintf("%.2f", myTicks),cex.axis=2,las=1)
  mtext(side=2,"Relative Intensity",col="black",cex = 2.5,line = 5.2)

  if (supress == FALSE) {
    x_range <- xlim[2] - xlim[1]
    y_position <- vector("numeric")
    if (num_identifications == 1) {
      y_position[i] <- aaa$expt_intensity[i]
    }
    else {
      for (i in 1:(num_identifications - 1)) {
        if ((aaa$expt_mz[i + 1] - aaa$expt_mz[i])/x_range <
            0.025 & all.equal(aaa$expt_intensity[i],
                              100) != TRUE) {
          y_position[i] <- aaa$expt_intensity[i] +
            40
          lines(rep(aaa$expt_mz[i], 2), c(aaa$expt_intensity[i],
                                          aaa$expt_intensity[i] + 80),
                lty = 5, col = color[i],lwd=2)
        }
        else y_position[i] <- aaa$expt_intensity[i]
      }
    }
    y_position[num_identifications] <- aaa$expt_intensity[num_identifications]
    text(aaa$expt_mz, y_position + 90, labels = aaa$ms2type,
         col = color, srt = 0, cex = 1.5)
  }
  seq_vector <- strsplit(as.character(aaa$ms1seq)[1],
                         split = "")[[1]]
  num_residues <- length(seq_vector)
  plot.window(xlim = c(1, 20), ylim = c(0, 10))
  text(c(1:num_residues), 9, labels = seq_vector,cex=2.5)
  lines(c(location[i] + 0.55, location[i] + 0.55),
        c(8.5,9.5), col = "black")

  for (i in 1:length(series)) {

    if (series[i] == "b" | series[i] == "c"){
      lines(c(location[i] + 0.55, location[i] + 0.55),
        c(8.5,9.5), col = "black")

      lines(c(location[i] + 0.55, location[i] + 0.95),
            c(9.5,9.5), col = "black")

    } else{lines(c(num_residues - location[i] + 0.55, num_residues -
                   location[i] + 0.55), c(8.5, 9.5), col = "black")

    lines(c(num_residues -location[i] + 0.55,  num_residues -location[i] + 0.15),
          c(8.5,8.5), col = "black")}

  }
  text(18, 9, label)

  }}



  }
  return(aaa)
}
