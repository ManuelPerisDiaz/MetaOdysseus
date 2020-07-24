#' Calculates permuted protein sequences according to number of Cys residues present.
#'
#' @param file peak file
#' @param do.plot do plot
#' @param t1 mz left
#' @param t2 mz right
#' @param dir directory to save
#' @param csv save in csv fileformat
#' @param orig Used original file
#' @param nameTry name to save
#' @examples
#' @return Plot peak
#' @export
#'
plotPeakBot<-function (file, do.plot = TRUE, t1 = 0, t2 = 0,csv=FALSE,orig=FALSE,nameTry=nameTry,dir="C:")
{
  oldpar <- par(no.readonly =TRUE)
  on.exit(par(oldpar))
  if (!is.logical(do.plot)) {
    stop("do.plot parameter must be logical")
  }
  if (t1 < 0 | t2 < 0 | t1 > t2) {
    stop("No possible values for t1 and t2")
  }


  if(csv== "TRUE"){

    peak <- read.csv(file, header = FALSE)}

  if(orig=="FALSE"){peak <-   file[,1]}else{
    peak <-  file$OR

  }


  if (ncol(peak) != 2) {
    stop("The file must have two columns: x and y")
  }
  if (t2 > t1) {
    peak <- peak[which(peak[, 1] > t1 & peak[, 1] < t2),
                 ]
  }

  y <- peak[, 2]
  y = (y-min(y))/(max(y)-min(y))
  x <- peak[, 1]
  if (do.plot == TRUE) {
    opar <- par(no.readonly =TRUE)
      d <- ggplot(peak, aes(x = x, y = y)) + geom_line(color = "red",
                                                     size = 0.5, show.legend = FALSE)
    d2 <- d + theme_bw()
    d22 <- d2  + theme(plot.title = element_text(size = 30),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border    = element_rect(colour = "black"),
                       axis.title = element_text(size = 20), axis.title.y = element_text(margin = margin(0, 20, 0, 0)),axis.text.y = element_text(size=15,colour = "black"),axis.text.x = element_text(size=15,colour = "black"), legend.text = element_blank(),
                       legend.title = element_text(size = 15)) + scale_y_continuous(name = "Relative Intensity")


    d22 <- d22 + xlab("m/z")+ggtitle(paste(as.character(nameTry)))
    on.exit(par(opar))
     ggsave(d22,filename=file.path(dir,paste("PMF_",nameTry,".tiff",sep="")), dpi = 800, width = 6, height = 4, units = 'in')

    }
  return(peak)
}


