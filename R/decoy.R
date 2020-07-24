#' Produce a decoy database
#'
#' @param sequence protein sequence
#' @param Npermu Number of permutations
#' @param df_full target database
#' @examples
#' @return Decoy database
#' @export

decoy<-function(sequence,Npermu=1,df_full){
  peptideVector <- strsplit(sequence, split = "")[[1]]
  Npermu=Npermu

  Permuted_fasta_lista<-c()
  for(i in 1:Npermu){
    Permuted_fasta <- permute(peptideVector)
    Permuted_fasta_lista[i]<-list(Permuted_fasta)
  }


  for(i in 1:length(Permuted_fasta_lista)){
    sequence<-paste(Permuted_fasta_lista[[i]],collapse="")
    decoy<-msms.batch(sequence,mzML = TRUE,missed=2,label="decoy",tol=2, Score=1)
  }

  decoy<-Filter(function(x) all(!is.na(x)),decoy)

#  df <- data.frame(matrix(unlist(decoy), nrow=length(decoy), byrow=T))
  df <- data.frame(matrix(unlist(data.frame(decoy)), nrow=length(decoy), byrow=T), stringsAsFactors=FALSE)
  colnames(df) <- c("Sequence","Mean Error (Da)","N. ions matched","Per. peaks matched","Per. ions matched","Int. ions matched","Score","Parent ion")

  decoy_full<-cbind(df,label="decoy")



  full_target_decoy<-rbind(decoy_full,df_full)

  colnames(full_target_decoy) <- c("Sequence","Mean Error (Da)","N. ions matched","Per. peaks matched","Per. ions matched","Int. ions matched","Score","Parent ion","label")

  full_target_decoy<-full_target_decoy[order(as.numeric(full_target_decoy$Score)),]

  write.csv2(full_target_decoy,file="CysmPro_annotated_target_decoy_nofilter.csv")

  return(full_target_decoy)

}
