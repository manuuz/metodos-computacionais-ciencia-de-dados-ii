pdf_cauchy <- function(x){
  result <- 1 / (pi * (1 + x^2))
  return (result)
}

mh_cauchy <- function(n, x0, sigma){
  chain <- numeric(n)
  chain[1] <- x0
  accepted = 0
  
  for (i in 2:length(chain)){
    x_current <- chain[i-1]
    
    # proposta: Normal com média de valor atual e sigma de entrada
    y <- rnorm(1, x_current, sigma)
    
    # R = f(Y) / f(x(t)) -> a g(.) se anula por ser uma normal
    ratio <- pdf_cauchy(y) / pdf_cauchy(x_current)
    
    if (runif(1) < ratio){
      chain[i] <- y 
      accepted<- accepted + 1
    } else {
      chain[i] <- x_current
    }
  }
  
  ARate <- accepted / (n-1)
  
  return(list(chain = chain, ARate = ARate))
}

n <- 10000
burn_in <- 1000

set.seed(123)  
result <- mh_cauchy(n, x0 = 5, sigma = 2)


chain <- result$chain
post_burn_chain <- chain[(burn_in + 1):n]

cat("Taxa de aceitação:", round(result$ARate, 4), "\n")

plot(post_burn_chain, type = "l", col = "blue", 
     main = "Traço da Cadeia (após burn-in)",
     xlab = "Iteração", ylab = "Valor")

hist(post_burn_chain, breaks = 50, freq = FALSE, 
     main = "Histograma vs Densidade Teórica",
     xlab = "x", ylab = "Densidade", col = "lightblue")
curve(dcauchy(x, location = 0, scale = 1), add = TRUE, 
      col = "red", lwd = 2, lty = 2)

