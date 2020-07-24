#' Calculates permuted protein sequences according to number of Cys residues present.
#'
#' @param file peak file
#' @param do.plotB do plot
#' @param do.plotA do plot
#' @param t1 mz left
#' @param t2 mz right
#' @param axis name axis
#' @param norm Normalization intensity
#' @param ID Name of file saved
#' @param csv save in csv fileformat
#' @param dir Directory to save file
#' @examples
#' @return Plot peak
#' @export

plot_Peak<-function (file, do.plotB = TRUE,do.plotA = TRUE,axis=c("m/z"),norm=TRUE, t1 = 1, t2 = 10000,ID,csv=FALSE,dir=getwd())
{
  oldpar <- par(no.readonly =TRUE)
  on.exit(par(oldpar))

  if (t1 < 0 | t2 < 0 | t1 > t2) {
    stop("No possible values for t1 and t2")
  }


  if(csv== "TRUE"){

  peak <- read.csv(file, header = FALSE)}

peak <-   file




  peakB <- peak[which(peak[, 1] > t1 & peak[, 1] < t2),]

  y <- peakB[, 2]

  x <- peakB[, 1]


  if (norm=="TRUE"){
    y = (y-min(y))/(max(y)-min(y))
  }
  if(class(peakB)=="matrix"){
    peakB<-as.data.frame(peakB)

  }

  if (do.plotA == TRUE) {
    dev.off()
     opar <- par(no.readonly =TRUE)
    d <- ggplot(data.frame(peakB), aes(x = x, y = y)) + geom_line(color = "black",
                                                                 size = 0.4, show.legend = FALSE)
    d2 <- d + theme_bw()
    d22 <- d2  + theme(plot.title = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border    = element_rect(colour = "black"),
                       axis.title = element_text(size = 20), axis.title.y = element_text(margin = margin(0, 20, 0, 0)),axis.text.y = element_text(size=15,colour = "black"),axis.text.x = element_text(size=15,colour = "black"), legend.text = element_blank(),
                       legend.title = element_text(size = 15)) + scale_y_continuous(name = "Relative Abundance")

    d22 <- d22 + xlab(axis)
    on.exit(par(opar))

   ggsave(d22,filename=file.path(dir,paste("MSA", ID, ".tiff",sep="")), dpi = 800, width = 4, height = 4, units = 'in')

  }

  if (do.plotB == TRUE) {
    dev.off()
    opar <- par(no.readonly =TRUE)
    d <- ggplot(data.frame(peakB), aes(x = x, y = y)) + geom_line(color = "black",
                                                                  size = 0.4, show.legend = FALSE)
    d2 <- d + theme_classic()
    d22 <- d2  + theme(plot.title = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       axis.line=element_line(colour="black"),
                       axis.title = element_text(size = 20), axis.title.y = element_text(margin = margin(0, 20, 0, 0)),axis.text.y = element_text(size=15,colour = "black"),axis.text.x = element_text(size=15,colour = "black"), legend.text = element_blank(),
                       panel.border = element_blank(),  legend.title = element_text(size = 15)) + scale_y_continuous(name = "Relative Abundance")

    d22 <- d22 + xlab(axis)
    on.exit(par(opar))

    ggsave(d22,filename=file.path(dir,paste("MSB", ID, ".tiff",sep="")), dpi = 800, width = 12, height = 4, units = 'in')

  }

  return(peakB)
}



