testez <- function(x, mu0, sigma, m = 10000) {
  n <- length(x)
  ep_mu <- sigma / sqrt(n)
  z_obs <- (mean(x) - mu0) / ep_mu
  xbarra_estrela <- rnorm(m, mean = mu0, sd = ep_mu)
  
  # estatística para cada média simulada
  z_estrela <- (xbarra_estrela - mu0) / ep_mu
  pvalor_mc <- sum(abs(z_estrela) >= abs(z_obs)) / m
  
  return(list(
    z_obs = z_obs,
    pvalor_mc = pvalor_mc
  ))
}

library(TeachingDemos) # função z.test
n <- 30
sigma <- 5
mu0 <- 10
mu1 <- 11.5
alpha <- 0.05
m <- 10000

rejeicoes_mc <- numeric(m)
rejeicoes_norm <- numeric(m)

set.seed(2023032088)
for (i in 1:m) {
  x_amostra <- rnorm(n, mean = mu1, sd = sigma)
  
  resultado_mc <- testez(x_amostra, mu0 = mu0, sigma = sigma, m = m)
  if (resultado_mc$pvalor_mc < alpha) {
    rejeicoes_mc[i] <- 1
  }

  resultado_norm <- z.test(x_amostra, mu = mu0, stdev = sigma,
                           alternative = "two.sided")
  if (resultado_norm$p.value < alpha) {
    rejeicoes_norm[i] <- 1
  }
}

poder_mc <- mean(rejeicoes_mc)
poder_norm <- mean(rejeicoes_norm)

ep_poder_mc <- sqrt(poder_mc * (1 - poder_mc) / m)
ep_poder_norm <- sqrt(poder_norm * (1 - poder_norm) / m)

tabela <- data.frame(
  Teste = c("THMC (Monte Carlo)", "Normal Analítico (Z-Test)"),
  Poder_Estimado = c(poder_mc, poder_norm),
  Erro_Padrao_MC = c(ep_poder_mc, ep_poder_norm)
)

kable(tabela, caption = "Comparação do Poder de Teste (mu real = 11.5)")