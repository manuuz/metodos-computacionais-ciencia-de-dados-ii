
rweib <- function(n, alpha, beta){
  U <- runif(n)
  x <- (-log(U) / alpha)^(1/beta)
  return(x)
}

densidade_teorica <- function(x, alpha, beta) {
  alpha * beta * x^(beta - 1) * exp(-alpha * x^beta)
}

set.seed(123)  

alpha <- 2
beta <- 1.5
n <- 10000

amostras <- rweib(n, alpha, beta)

hist(amostras, breaks = 30, probability = TRUE, 
     main = "Histograma vs Densidade Teórica",
     xlab = "x", col = "lightblue")

x_seq <- seq(0, max(amostras), length.out = 1000)
lines(x_seq, densidade_teorica(x_seq, alpha, beta), col = "red", lwd = 2)

legend("topright", legend = c("Amostras", "Teórica"), 
       col = c("lightblue", "red"), lwd = c(2, 2))

