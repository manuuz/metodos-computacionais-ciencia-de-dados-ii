set.seed(2023032088)
m <- 5000
alpha <- 0.05
mu0 <- 500
sigma <- 100
n <- c(10, 20, 30, 40, 50)

mu_alternativa <- seq(400, 600, by = 5)

# matriz para armazenar as estimativas de poder (linhas = mu, colunas = n)
poder <- matrix(NA, nrow = length(mu_alternativa), ncol = length(n))
colnames(poder) <- paste0("n=", n)

for (j in 1:length(n)) {
  novo_n <- n[j]
  
  for (i in 1:length(mu_alternativa)) {
    mu_real <- mu_alternativa[i]
    
    pvalores <- replicate(m, expr = {
      x <- rnorm(novo_n, mean = mu_real, sd = sigma) 
      teste_t <- t.test(x, alternative = "two.sided", mu = mu0)
      teste_t$p.value
      }
    )
    poder[i, j] <- mean(pvalores < alpha)
  }
}

cores <- c("blue", "red", "green4", "purple", "orange")
curvas <- paste("n =", n)
par(mar = c(5, 4, 4, 6))

plot(mu_alternativa, poder[, 1], 
     type = "l", 
     col = cores[5], 
     lwd = 2,
     ylim = c(0, 1),
     xlab = expression(paste("Média Populacional (", mu, ")")), 
     ylab = "Poder",
     main = "Curva de Poder Empírica para o Teste t Bilateral",
     las = 1)

for (j in 2:length(n)) {
  lines(mu_alternativa, poder[, j], col = cores[j], lwd = 2)
}

# referência
abline(h = alpha, lty = 2, col = "gray50")
abline(v = mu0, lty = 3, col = "gray50")

legend("topright",
       inset = c(-0.3,0),
       legend = curvas,
       col = cores, 
       lwd = 2,
       bty = "n",
       xpd = TRUE)

kable(round(poder[c(1, 11, 21, 31, 41), ], 4), caption = "Poder estimado")