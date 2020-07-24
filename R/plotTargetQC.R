plotTargetQC<-function(FiltretedINt){
d <- ggplot(data.frame(FiltretedINt), aes(x = FiltretedINt$mz, y = FiltretedINt$Inten)) + geom_line(color = "red",
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
