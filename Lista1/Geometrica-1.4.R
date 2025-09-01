# Número de falhas antes do primeiro sucesso
rgeom <- function(n,p){
  u <- runif(n, 0, 1)
  q <- 1 - p
  
  k <- floor(log(u) / log(q))
  return(k)
}

n <- 10000
p <- 0.25

x <- rgeom(n, p)
hist(x,
     prob = TRUE,
     main = "Distribuição Geométrica (k >= 1)",
     xlab = "Número de Tentativas (k)",
     ylab = "Frequência Relativa")

x_vals <- 1:max(x)
y_vals <- dgeom(x_vals - 1, prob = p)
lines(x_vals, y_vals, col = "tomato3", lwd = 2)