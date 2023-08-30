#' @export

CondListToString <- function(namesList)
  return(lapply(namesList, function(x)  paste(x, collapse = ", ")))

GetRownamesCML <- function(firsts, seconds)
  return(mapply(function(x, y) str_c("[", x, "] vs. [", y, "]"), CondListToString(firsts), CondListToString(seconds)))

BonferroniAndFilter <- function(markers, nTests){
  #Correcting multiple condition testing with Bonferroni
  markers$p_val_adj <- markers$p_val_adj * nTests
  markers <- subset(markers, p_val_adj < 0.05)
  return(markers)
}

GetA13ANames <- function(){
  firstsA13A <- list("Activin A", "Activin A and I-BRD9", "Activin A", "SB-431542", "Activin A and I-BRD9", "SB-431542",
                     "Activin A", c("Activin A and I-BRD9", "SB-431542"), "Activin A and I-BRD9", c("Activin A", "SB-431542"), "SB-431542",
                     c("Activin A", "Activin A and I-BRD9"))
  secondsA13A <- firstsA13A[c(2,1,4,3,6,5,8,7,10,9,12,11)]
  RownamesCMLA13A <- GetRownamesCML(firstsA13A, secondsA13A)
  return (list(firstsA13A, secondsA13A, RownamesCMLA13A))
}

GetPatientNames <- function(){
  firstsPatient <- list("DMSO", "I-BRD9", "DMSO", "Gemcitabine", "DMSO", "I-BRD9 and gemcitabine", "I-BRD9",
                        "Gemcitabine", "I-BRD9", "I-BRD9 and gemcitabine", "Gemcitabine", "I-BRD9 and gemcitabine",
                        c("DMSO", "Gemcitabine"), c("I-BRD9", "I-BRD9 and gemcitabine"))
  secondsPatient <- firstsPatient[c(2,1,4,3,6,5,8,7,10,9,12,11,14,13)]
  RownamesCMLPatient <- GetRownamesCML(firstsPatient, secondsPatient)
  return (list(firstsPatient, secondsPatient, RownamesCMLPatient))
}




