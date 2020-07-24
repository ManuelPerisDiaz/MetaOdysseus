#' Targeted MS/MS annotation
#'
#' @param sequence Protein sequence fasta
#' @param fragments fragments
#' @param IAA IAA
#' @param N15 N15
#' @param custom custom
#' @examples
#' @return Targeted MS/MS annotation
#' @export
MSMS<-function (sequence, fragments = "by", IAA = TRUE, N15 = FALSE,
          custom = list())
{

    results_list <- vector("list")
    for (sequence_number in 1:length(sequence)) {
    peptide_vector <- strsplit(sequence[sequence_number],
                               split = "")[[1]]
    peptide_length <- length(peptide_vector)
    if (peptide_length < 2)
      stop("sequence must contain two or more residues")
    C <- 12
    H <- 1.0078250321
    O <- 15.9949146221
    S <- 31.97207069
    N <- ifelse(N15 == TRUE, 15.0001088984, 14.0030740052)
    proton <- 1.007276466
    electron <- 0.00054857990943

    residueMass <- function(residue) {
      if (residue == "A")
        mass = C * 3 + H * 5 + N + O
      if (residue == "R")
        mass = C * 6 + H * 12 + N * 4 + O
      if (residue == "N")
        mass = C * 4 + H * 6 + N * 2 + O * 2
      if (residue == "D")
        mass = C * 4 + H * 5 + N + O * 3
      if (residue == "E")
        mass = C * 5 + H * 7 + N + O * 3
      if (residue == "Q")
        mass = C * 5 + H * 8 + N * 2 + O * 2
      if (residue == "G")
        mass = C * 2 + H * 3 + N + O
      if (residue == "H")
        mass = C * 6 + H * 7 + N * 3 + O
      if (residue == "I")
        mass = C * 6 + H * 11 + N + O
      if (residue == "L")
        mass = C * 6 + H * 11 + N + O
      if (residue == "K")
        mass = C * 6 + H * 12 + N * 2 + O
      if (residue == "M")
        mass = C * 5 + H * 9 + N + O + S
      if (residue == "F")
        mass = C * 9 + H * 9 + N + O
      if (residue == "P")
        mass = C * 5 + H * 7 + N + O
      if (residue == "S")
        mass = C * 3 + H * 5 + N + O * 2
      if (residue == "T")
        mass = C * 4 + H * 7 + N + O * 2
      if (residue == "W")
        mass = C * 11 + H * 10 + N * 2 + O
      if (residue == "Y")
        mass = C * 9 + H * 9 + N + O * 2
      if (residue == "V")
        mass = C * 5 + H * 9 + N + O
      if (residue == "C" )
        mass = C * 3 + H * 5 + N + O + S
      if (residue == "X" & IAA == TRUE)
        mass <- ifelse(N15 == FALSE, C * 5 + H * 8 +
                         N * 2 + O * 2 + S, C * 5 + H * 8 + N + 14.0030740052 +
                         O * 2 + S)
      if (length(custom) != 0)
        for (i in 1:length(custom$code)) if (residue ==
                                             custom$code[i])
          mass = custom$mass[i]
        return(mass)
    }
    masses <- sapply(peptide_vector, residueMass)
    pm <- sum(masses)
    p1 <- round(pm + H * 2 + O + proton, digits = 3)
    p2 <- round((pm + H * 2 + O + (2 * proton))/2, digits = 3)
    p3 <- round((pm + H * 2 + O + (3 * proton))/3, digits = 3)
    if (fragments == "by") {
      b1 <- vector(mode = "numeric", length = 0)
      b2 <- vector(mode = "numeric", length = 0)
      bs <- vector(mode = "character", length = 0)
      bi <- vector(mode = "integer", length = 0)
      y1 <- vector(mode = "numeric", length = 0)
      y2 <- vector(mode = "numeric", length = 0)
      ys <- vector(mode = "character", length = 0)
      yi <- vector(mode = "integer", length = 0)
      for (i in 1:(peptide_length - 1)) {
        mass <- sum(masses[1:i])
        b1[i] <- round(mass + proton, digits = 3)
        b2[i] <- round((b1[i] + proton)/2, digits = 3)
        bs[i] <- paste(peptide_vector[1:i], collapse = "")
        bi[i] <- i
      }
      for (j in 2:peptide_length) {
        mass <- sum(masses[j:peptide_length])
        y1[j - 1] <- round(mass + H * 2 + O + proton,
                           digits = 3)
        y2[j - 1] <- round((y1[j - 1] + proton)/2, digits = 3)
        ys[j - 1] <- paste(peptide_vector[j:peptide_length],
                           collapse = "")
        yi[j - 1] <- peptide_length - j + 1
      }
      ms1seq <- rep(sequence[sequence_number], times = ((2 *
                                                           (length(bi))) + (2 * (length(yi)))))
      ms1z1 <- rep(p1, times = ((2 * (length(bi))) + (2 *
                                                        (length(yi)))))
      ms1z2 <- rep(p2, times = ((2 * (length(bi))) + (2 *
                                                        (length(yi)))))
      ms1z3 <- rep(p3, times = ((2 * (length(bi))) + (2 *
                                                        (length(yi)))))
      ms2seq <- c(rep(bs, times = 2), rep(ys, times = 2))
      b1.type <- paste("[b", bi, "]1+", sep = "")
      b2.type <- paste("[b", bi, "]2+", sep = "")
      y1.type <- paste("[y", yi, "]1+", sep = "")
      y2.type <- paste("[y", yi, "]2+", sep = "")
      ms2type <- c(b1.type, b2.type, y1.type, y2.type)
      ms2mz <- c(b1, b2, y1, y2)
    }
    if (fragments == "cz") {
      c1 <- vector(mode = "numeric", length = 0)
      c2 <- vector(mode = "numeric", length = 0)
      cs <- vector(mode = "character", length = 0)
      ci <- vector(mode = "integer", length = 0)
      z1 <- vector(mode = "numeric", length = 0)
      z2 <- vector(mode = "numeric", length = 0)
      zs <- vector(mode = "character", length = 0)
      zi <- vector(mode = "integer", length = 0)
      for (i in 1:(peptide_length - 1)) {
        mass <- sum(masses[1:i])
        c1[i] <- round(mass + 3 * H + N + proton, digits = 3)
        c2[i] <- round((c1[i] + proton)/2, digits = 3)
        cs[i] <- paste(peptide_vector[1:i], collapse = "")
        ci[i] <- i
      }
      for (j in 2:peptide_length) {
        mass <- sum(masses[j:peptide_length])
        z1[j - 1] <- round(mass + O - N, digits = 3)
        z2[j - 1] <- round((z1[j - 1] + proton)/2, digits = 3)
        zs[j - 1] <- paste(peptide_vector[j:peptide_length],
                           collapse = "")
        zi[j - 1] <- peptide_length - j + 1
      }
      ms1seq <- rep(sequence[sequence_number], times = ((2 *
                                                           (length(ci))) + (2 * (length(zi)))))
      ms1z1 <- rep(p1, times = ((2 * (length(ci))) + (2 *
                                                        (length(zi)))))
      ms1z2 <- rep(p2, times = ((2 * (length(ci))) + (2 *
                                                        (length(zi)))))
      ms1z3 <- rep(p3, times = ((2 * (length(ci))) + (2 *
                                                        (length(zi)))))
      ms2seq <- c(rep(cs, times = 2), rep(zs, times = 2))
      c1.type <- paste("[c", ci, "]1+", sep = "")
      c2.type <- paste("[c", ci, "]2+", sep = "")
      z1.type <- paste("[z", zi, "]1+", sep = "")
      z2.type <- paste("[z", zi, "]2+", sep = "")
      ms2type <- c(c1.type, c2.type, z1.type, z2.type)
      ms2mz <- c(c1, c2, z1, z2)
    }


    results_list[[sequence_number]] <- data.frame(ms1seq,
                                                  ms1z1, ms1z2, ms1z3, ms2seq, ms2type, ms2mz)
    }

      return(as.data.frame(do.call("rbind", results_list)))
}
