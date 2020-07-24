#' Permuted protein sequence.
#'
#' Calculates permuted protein sequences according to number of Cys residues present.
#'
#' @param sequence Protein sequence.
#' @param NModif Number of total modifications expected in the Cys-protein residues.
#' @examples
#' \dontrun{
#' sequence<-c("KSCCSCCPAECEK")
#' listSeq<-permuSEQ(sequence,NModif=4)
#' }
#' @return Permuted Cys-protein sequences.
#' @export

permuSEQ<-function(sequence,NModif=3){
  results_list <- vector("list")
  grepEx<-which(strsplit(sequence,"")[[1]]=="C")
  peptide_vector <- strsplit(sequence, split = "")[[1]]
  peptide_length <- length(peptide_vector)
  NCys<-sum(as.numeric(peptide_vector=="C"))
  comb<-combinations(n=NCys,r=NModif,v=grepEx,set=T,repeats.allowed=F)

  M2<-c()
  M<-c()
  for (i in 1:nrow(comb)){
    peptide_vector <- strsplit(sequence,
                               split = "")[[1]]
    peptide_vector[comb[i,]]<-"X"
    M<-peptide_vector
    M2[i]<-list(M)
  }

  f<-c()

  for(i in 1:length(M2)){
    result<-c()
    result<-paste(M2[[i]], collapse="")
    f[i]<-list(result)
  }

  return(f)
}
