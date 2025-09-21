# estruturação
library(knitr)
set.seed(2023032088)
m <- 1000
n <- c(10, 30, 100)
resultados <- list()

simulacao <- function(n, m, dist, mu, sigma2, r_function) {
  sigma2_real <- sigma2
  nao_viciado <- numeric(m)
  viciado <- numeric(m)
  
  for (j in 1:m) {
    x <- r_function(n)
    # cálculo dos estimadores
    nao_viciado[j] <- var(x) 
    viciado[j] <- (n - 1) / n * nao_viciado[j]
  }
  
  vies <- function(estimador) {
    G_vies <- estimador - sigma2_real 
    vies_estimado <- mean(G_vies)
    ep_vies <- sd(G_vies) / sqrt(m)
    return(c(Vies = vies_estimado, SE_Vies = ep_vies))
  }
  
  eqm <- function(estimador) {
    G_eqm <- (estimador - sigma2_real)^2 
    eqm_estimado <- mean(G_eqm)
    ep_eqm <- sd(G_eqm) / sqrt(m)
    return(c(EQM = eqm_estimado, SE_EQM = ep_eqm))
  }

  vies1 <- vies(nao_viciado)
  eqm1  <- eqm(nao_viciado)
  vies2 <- vies(viciado)
  eqm2  <- eqm(viciado)
  
  res1 <- data.frame(
    Distribuicao = dist, n = n, Estimador = "Sigma1^2 (Não Viciado)",
    vies = vies1["Vies"], se_vies = vies1["SE_Vies"],
    eqm = eqm1["EQM"], se_eqm = eqm1["SE_EQM"]
  )
  res2 <- data.frame(
    Distribuicao = dist, n = n, Estimador = "Sigma2^2 (Viciado)",
    vies = vies2["Vies"], se_vies = vies2["SE_Vies"],
    eqm = eqm2["EQM"], se_eqm = eqm2["SE_EQM"]
  )
  
  return(rbind(res1, res2))
}

# exemplos
# N(10, 4), sigma^2 = 4
mu_norm <- 10
sigma2_norm <- 4
r_norm <- function(n) { rnorm(n, mean = mu_norm, sd = sqrt(sigma2_norm)) }

for (n in n) {
  resultados[[length(resultados) + 1]] <- simulacao(
    n = n, m = m, dist = "Normal (N(10, 4))", mu = mu_norm, 
    sigma2 = sigma2_norm, r_function = r_norm
  )
}

# Exp(0.5), sigma^2 = 1/0.5^2 = 4
lambda_exp <- 0.5
mu_exp <- 1/lambda_exp
sigma2_exp <- 1/lambda_exp^2
r_exp <- function(n) { rexp(n, rate = lambda_exp) }

for (n in n) {
  resultados[[length(resultados) + 1]] <- simulacao(
    n = n, m = m, dist = "Exponencial (Exp(0.5))", mu = mu_exp, 
    sigma2 = sigma2_exp, r_function = r_exp
  )
}


tabela_final <- do.call(rbind, resultados)
tabela_final[, 4:7] <- round(tabela_final[, 4:7], 6)
kable(tabela_final, caption = "Viés e EQM da Variância")