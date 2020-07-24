


target.MSMS.Multiple<-function(theoryMultiple,Deconv_CysMpro,permutation=FALSE,native=FALSE,t=0.4,b=0,plot=FALSE,supress=FALSE,xlim = c(100, 8000)){

  final<-c()
  for (i in 1:length(theoryMultiple)){
    results<-c()

    if(native==TRUE){

      for (z in 1:(length((theoryMultiple)[[i]]))){

        if(plot==TRUE){
        tiff(paste0("theoryMultiple",z,".tiff",sep=""),width = 10, height = 6, units = 'in',res=600)
        }

        if(permutation==TRUE){
          results<-targetMSMS_permu(Deconv_CysMpro,theoryMultiple[[i]][z],t=t,b=b,xlim =xlim,plot=plot,supress=supress)


        }else{

          results<-targetMSMS(Deconv_CysMpro,theoryMultiple[[i]][z],mzw=t,b=b,xlim =xlim,plot=plot,supress=supress)

        }


      if(plot==TRUE){
       dev.off()}
       final<-rbind(final,results)


      }

    if(!is.null(final)){
         write.csv2(final,file=paste("CID_",i,".csv",sep=""))
      }

    }else{

      tiff(paste0("theoryMultiple",i,".tiff",sep=""),width = 10, height = 6, units = 'in',res=600)
     results<-targetMSMS(Deconv_CysMpro,theoryMultiple[[i]],mzw=t,b=b,xlim =xlim)
  dev.off()
  final[i]<-list(results)
  write.csv2(final[i],paste0("Sequence",i,".csv",sep=""))

   }}

  return(final)
}


