plot_processed<-function (file,norm=FALSE, t1 = 0, t2 = 0,csv=FALSE,dir=getwd())
{

  oldpar <- par(no.readonly =TRUE)
  on.exit(par(oldpar))
  if (t1 < 0 | t2 < 0 | t1 > t2) {
    stop("No possible values for t1 and t2")
  }


  if(csv== "TRUE"){

    peak <- read.csv(file, header = FALSE)}else{

      peak<-file
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
    y = (y-min(y))/(max(y)-min(y))
  }

  if(class(peak)=="matrix"){
    peakB<-as.data.frame(peak)

  }

    opar <- par(no.readonly =TRUE)
    d <- ggplot(data.frame(peak), aes(x = x, y = y)) + geom_line(color = "red",
                                                                 size = 0.5, show.legend = FALSE)
    d2 <- d + theme_bw()
    d22 <- d2  + theme(plot.title = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border    = element_rect(colour = "black"),
                       axis.title = element_text(size = 20), axis.title.y = element_text(margin = margin(0, 20, 0, 0)),axis.text.y = element_text(size=15,colour = "black"),axis.text.x = element_text(size=15,colour = "black"), legend.text = element_blank(),
                       legend.title = element_text(size = 15)) + scale_y_continuous(name = "Relative Intensity")

    d22 <- d22 + xlab("m/z")
    on.exit(par(opar))
    ggsave(d22,filename=file.path(dir,"MS_processed.tiff"), dpi = 800, width = 6, height = 4, units = 'in')



  return(peak)

}
