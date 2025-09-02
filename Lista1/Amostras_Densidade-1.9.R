rfdp <- function(n){
  U <- runif(n)
  
  # Inversa da FDA F^-1(u)
  amostras <- 4 * sqrt(U)
  
  return(amostras)
}

n <- 10000

amostra_gerada <- rfdp(n)

hist(amostra_gerada, 
     breaks = 30, 
     freq = FALSE, 
     main = "Amostra Gerada vs. Densidade Teórica",
     xlab = "x",
     ylab = "Densidade",
     col = "lightblue",
     border = "white")

curve(x/8, 
      from = 0, 
      to = 4, 
      add = TRUE, 
      col = "tomato3", 
      lwd = 2) 
