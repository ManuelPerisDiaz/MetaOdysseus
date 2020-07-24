plotESI<-function (file,norm=FALSE, do.plot = TRUE, t1 = 0, t2 = 0,csv=FALSE,orig=TRUE,QC1=TRUE,QC2=TRUE)
{

  if (!is.logical(do.plot)) {
    stop("do.plot parameter must be logical")
  }
  if (t1 < 0 | t2 < 0 | t1 > t2) {
    stop("No possible values for t1 and t2")
  }


  if(csv== "TRUE"){

    peak <- read.csv(file, header = FALSE)}

  if(orig=="FALSE"){peak <-   file}else{
    peak <-  file$OR

  }


  if (t2 > t1) {
    peak <- peak[which(peak[, 1] > t1 & peak[, 1] < t2),
                  ]
    y <- peak[, 2]

    x <- peak[, 1]

    }

  y <- peak[, 2]

  x <- peak[, 1]

  if (norm=="TRUE"){
  y = (y-min(y))/(max(y)-min(y))}

  if(class(peak)=="matrix"){
    peakB<-as.data.frame(peak)

  }

  if (do.plot == TRUE ) {

    d <- ggplot(data.frame(peak), aes(x = x, y = y)) + geom_line(color = "red",
                                                      size = 0.5, show.legend = FALSE)
    d2 <- d + theme_bw()
    d22 <- d2  + theme(plot.title = element_blank(),
                       axis.title = element_text(size = 30), axis.title.y = element_text(margin = margin(0, 20, 0, 0)),axis.text.y = element_text(size=15),axis.text.x = element_text(size=15), legend.text = element_blank(),
                       legend.title = element_text(size = 15)) + scale_y_continuous(name = "Relative Intensity")

    d22 <- d22 + xlab("m/z")
    d22 + theme_minimal()
    d22 + scale_x_continuous(expression(m/z))
    print(d22)
    ggsave(d22,filename="CysPro.tiff", dpi = 800, width = 6, height = 4, units = 'in')

  }

  if(QC1==TRUE){

    d <- ggplot(data.frame(peak), aes(x = x,xend=x,y=0, yend = y)) + geom_segment(color = "red",
                                                                 size = 0.5, show.legend = FALSE)
    d2 <- d + theme_bw()
    d22 <- d2  + theme(plot.title = element_blank(),
                       axis.title = element_text(size = 30), axis.title.y = element_text(margin = margin(0, 20, 0, 0)),axis.text.y = element_text(size=15),axis.text.x = element_text(size=15), legend.text = element_blank(),
                       legend.title = element_text(size = 15)) + scale_y_continuous(name = "Relative Intensity")

    d22 <- d22 + xlab("m/z")
    d22 + theme_minimal()
    d22 + scale_x_continuous(expression(m/z))
    print(d22)
    ggsave(d22,filename="QC1.tiff", dpi = 800, width = 6, height = 4, units = 'in')

  }

  if(QC2==TRUE){


    d <- ggplot(data.frame(peak), aes(x = x, y = y)) + geom_point(color = "red",
                                                                                  size = 2, show.legend = FALSE)
    d2 <- d + theme_bw()
    d22 <- d2  + theme(plot.title = element_blank(),
                       axis.title = element_text(size = 30), axis.title.y = element_text(margin = margin(0, 20, 0, 0)),axis.text.y = element_text(size=15),axis.text.x = element_text(size=15), legend.text = element_blank(),
                       legend.title = element_text(size = 15)) + scale_y_continuous(name = "Relative Intensity")

    d22 <- d22 + xlab("m/z deviation (Da)")
    d22 + theme_minimal()
    d22 + scale_x_continuous(expression(m/z))
    print(d22)
    ggsave(d22,filename="QC2.tiff", dpi = 800, width = 6, height = 4, units = 'in')

  }


  return(peak)


  }

