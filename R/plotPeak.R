#' Calculates permuted protein sequences according to number of Cys residues present.
#'
#' @param file peak file
#' @param do.plot do plot
#' @param t1 mz left
#' @param t2 mz right
#' @param csv save in csv fileformat
#' @param orig Used original file
#' @examples
#' @return Plot peak
#' @export
plotPeak<-function (file,do.plot = TRUE, t1 = 0, t2 = 0,csv=FALSE,orig=TRUE)
{

  if (!is.logical(do.plot)) {
    stop("do.plot parameter must be logical")
  }
  if (t1 < 0 | t2 < 0 | t1 > t2) {
    stop("No possible values for t1 and t2")
  }


  if(csv== "TRUE"){

  peak <- read.csv(file, header = FALSE)}

  if(orig=="FALSE"){peak <-   file$IN.Mz}else{
    peak <-  file$OR

  }


  if (t2 > t1) {
    peakB <- peak[which(peak[, 1] > t1 & peak[, 1] < t2),]
  }
  y <- peakB[, 2]

  x <- peakB[, 1]
  y = (y-min(y))/(max(y)-min(y))

  if(class(peakB)=="matrix"){
    peakB<-as.data.frame(peakB)

  }

  if (do.plot == TRUE) {

    d <- ggplot(peakB, aes(x = x, y = y)) + geom_line(color = "red",
                                                      size = 0.5, show.legend = FALSE)
    d2 <- d + theme_bw()
    d22 <- d2  + theme(plot.title = element_blank(),
                                          axis.title = element_text(size = 30), axis.title.y = element_text(margin = margin(0, 20, 0, 0)),axis.text.y = element_text(size=15),axis.text.x = element_text(size=15), legend.text = element_blank(),
                                          legend.title = element_text(size = 15)) + scale_y_continuous(name = "Relative intensity")

     d22 <- d22 + xlab("m/z")
    d22 + theme_minimal()
    d22 + scale_x_continuous(expression(m/z))
    print(d22)
   ggsave(d22,filename="CysPro.tiff", dpi = 800, width = 6, height = 4, units = 'in')

  }
  return(peakB)
}

