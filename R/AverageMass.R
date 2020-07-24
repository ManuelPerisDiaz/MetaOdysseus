AverageMass<-function (formula = list(), isotopes = list(), charge = 0) 
{
  defaultFormula <- list(C = 0, H = 0, N = 0, O = 0, S = 0, 
                         P = 0, Br = 0, Cl = 0, F = 0, Si = 0, x = 0)
  defaultFormula[names(formula)] <- formula
  defaultIsotopes <- list(C = 12.010736, H = 1.007941, N = 14.006703, 
                          O = 15.999405, S = 32.064787, P = 30.973762, 
                          Br = 79.903528, Cl = 35.452938, F = 18.998403, Si = 28.085499, 
                          x = 0)
  defaultIsotopes[names(isotopes)] <- isotopes
  if (charge < 0 & abs(charge) > defaultFormula$H) 
    stop("the number of negative charges exceeds the number of hydrogens in the formula list")
  mass <- (defaultFormula$C * defaultIsotopes$C + defaultFormula$H * 
             defaultIsotopes$H + defaultFormula$N * defaultIsotopes$N + 
             defaultFormula$O * defaultIsotopes$O + defaultFormula$S * 
             defaultIsotopes$S + defaultFormula$P * defaultIsotopes$P + 
             defaultFormula$Br * defaultIsotopes$Br + defaultFormula$Cl * 
             defaultIsotopes$Cl + defaultFormula$F * defaultIsotopes$F + 
             defaultFormula$Si * defaultIsotopes$Si + defaultFormula$x * 
             defaultIsotopes$x)
  if (charge != 0) 
    mass <- abs((mass + charge * 1.007276466)/charge)
  return(mass)
}