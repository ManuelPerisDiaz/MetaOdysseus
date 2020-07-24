PerPeptidefrag<-function(exp,target,b=5){
  result<-c()
  expt <- data.frame(mz = exp[, 1], intensity = exp[, 2])
  expt$normalized <- ((expt$intensity-min(expt$intensity))/(max(expt$intensity)-min(expt$intensity)))*100
  expt <- subset(expt, expt$normalized >= b) 
  result<-(sum(target$expt_intensity))/(sum(expt$normalized))*100
  return(result)
}
