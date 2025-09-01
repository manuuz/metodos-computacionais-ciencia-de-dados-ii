permutacao <- function(vetor){
  vetor_aux <- vetor
  n <- length(vetor_aux)
  
  for (k in n:2) {
    i <- floor(runif(1, 0, 1) * k) + 1

    temp <- vetor_aux[i]
    vetor_aux[i] <- vetor_aux[k]
    vetor_aux[k] <- temp
  }
  
  return(vetor_aux)
  
}
vetor <- 1:10
funcao <- permutacao(vetor)
funcao


