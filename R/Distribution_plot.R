#' Distribution plot
#'
#' @param pmfscore_totalb Output of the scoring
#' @param threshold threshold for scoring
#' @param binwidth Binwidth
#' @param sd sd
#' @param name name for saving
#' @param do.plotB If TRUE returns plot
#' @examples
#' @return Distribution plot for scored PMF MALDI-MS.
#' @export

Distribution_plot<-function(pmfscore_totalb,threshold,binwidth=0.2,
                            sd,name,do.plotB=FALSE){



d<-qplot(pmfscore_totalb,xlim=c((min(pmfscore_totalb)-1),max(threshold+sd+2)), binwidth=binwidth,geom="histogram",

      col=I ("black"), fill=I("white"),xlab = "Score")



d2<-d+    theme_bw() + theme(axis.title = element_text(size = 28), axis.text=element_text(size=24)
        ,axis.title.y=element_text(margin=margin(0,15,0,0)),

        axis.title.x=element_text(margin=margin(15,0,0,0)),

        legend.title = element_text(size = 22),legend.position="none")+scale_y_continuous(expression(Frequency))


d3<-d2+geom_segment(aes( x= threshold-sd,xend=threshold+sd,y=1,yend=1),size=1,  lineend = "square")+
  geom_point(aes(y=1,x=threshold))

d3

file_name = paste("Permutation_", name, ".tiff", sep="")
tiff(file_name, res = 800, width = 6, height = 4, units = 'in')
print(d3)
dev.off()

if (do.plotB == TRUE) {

  d<-qplot(pmfscore_totalb,xlim=c((min(pmfscore_totalb)-1),max(threshold+1)), binwidth=binwidth,geom="histogram",

           col=I ("black"), fill=I("white"),xlab = "Score")


  d2 <- d +   theme_void()+

    theme(axis.title = element_text(size = 28), axis.ticks.y=element_blank(),   axis.text.y=element_blank(), axis.text.x = element_text(size=24,color = "black"),
                               axis.title.y=element_text(margin=margin(0,15,0,0)),
                               panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                               axis.title.x=element_text(margin=margin(5,0,0,0)),
          panel.border = element_blank(),
                               legend.title = element_text(size = 22),legend.position="none")

  d2

  d3<-d2+geom_segment(aes( x= threshold-sd,xend=threshold+sd,y=1,yend=1),size=1,  lineend = "square")+
    geom_point(aes(y=1,x=threshold))

  d3

  file_name = paste("Permutation_", name, ".tiff", sep="")
  tiff(file_name, res = 800, width = 6, height = 4, units = 'in')
  print(d3)
  dev.off()


  }




}


