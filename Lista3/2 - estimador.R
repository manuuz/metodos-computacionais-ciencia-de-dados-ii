# parâmetros
K <- 10000 # amostras por replicação
M <- 200  # número de replicações
mu_Y <- 1
sigma2_Y <- 1
sigma2_XdadoY <- 4
theta_exato <- 0.5

estimativas_mc <- numeric(M)

set.seed(123456789)

for (m in 1:M) {
  Y <- rnorm(K, mean = mu_Y, sd = sqrt(sigma2_Y))
  epsilon <- rnorm(K, mean = 0, sd = sqrt(sigma2_XdadoY))
  
  # X condicional a Y
  X <- Y + epsilon
  
  estimativas_mc[m] <- mean(X > 1)
}

estimativas_ec <- numeric(M)

for (m in 1:M) {
  Y <- rnorm(K, mean = mu_Y, sd = sqrt(sigma2_Y))
  
  #condicional g(Y) = P(X>1|Y)
  g_Y <- pnorm((Y - 1) / 2)
  
  estimativas_ec[m] <- mean(g_Y)
}

estimativas_vc <- numeric(M)

for (m in 1:M) {
  Y <- rnorm(K, mean = mu_Y, sd = sqrt(sigma2_Y))
  
  # condicional e coeficiente ótimo
  g_Y <- pnorm((Y - 1) / 2)
  c_chapeu <- cov(g_Y, Y) / var(Y)
  
  Y_barra <- mean(Y)
  estimativas_vc[m] <- mean(g_Y) - c_chapeu * (Y_barra - mu_Y)
}

ep_mc <- sd(estimativas_mc)
ep_ec <- sd(estimativas_ec)
ep_vc <- sd(estimativas_vc)

media_mc <- mean(estimativas_mc)
media_ec <- mean(estimativas_ec)
media_vc <- mean(estimativas_vc)

# redução de variância (em comparação com o método simples)
reducao_ec <- (ep_mc^2 - ep_ec^2) / ep_mc^2 * 100
reducao_vc <- (ep_mc^2 - ep_vc^2) / ep_mc^2 * 100

resultados <- data.frame(
  `Método` = c("MC Simples", "Esperança Condicional", "EC + Variável de Controle"),
  `Média Estimada` = c(media_mc, media_ec, media_vc),
  `Erro Padrão Estimado` = c(ep_mc, ep_ec, ep_vc),
  `Redução de Variância` = c(0.0, reducao_ec, reducao_vc)
)
kable(resultados)

boxplot(estimativas_mc, estimativas_ec, estimativas_vc,
        names = c("MC Simples", "E. Condicional", "EC + VC"),
        main = "Variabilidade dos Estimadores (M=200)",
        ylab = "Estimativas de theta",
        col = c("lightblue", "lightgreen", "salmon"))
abline(h = theta_exato, col = "red", lty = 2, lwd = 2)