# ==============================================================================
# ESTUDO MC: VIÉS E EQM DOS ESTIMADORES DE VARIÂNCIA
# (Adaptado do Exercício 2.2)
# ==============================================================================

set.seed(42)

# 1. Parâmetros da Simulação
M <- 1000          # Número de réplicas Monte Carlo [2]
n_valores <- c(10, 30, 100)
resultados <- list()

# Função para realizar a simulação para uma dada distribuição
simulacao_variancia_mc <- function(n, M, dist_name, mu, sigma2, r_function) {
  
  # Variância populacional verdadeira
  sigma2_real <- sigma2
  
  # Matrizes para armazenar os estimadores
  sigma1_sq_hats <- numeric(M)
  sigma2_sq_hats <- numeric(M)
  
  for (j in 1:M) {
    # 1. Gerar a amostra
    x <- r_function(n)
    
    # 2. Calcular a variância amostral (n-1 no denominador)
    # R utiliza var(x) = sigma1_sq_hat
    sigma1_sq_hats[j] <- var(x) 
    
    # 3. Calcular o estimador viciado (n no denominador)
    # sigma2_sq_hat = (n-1)/n * sigma1_sq_hat
    sigma2_sq_hats[j] <- (n - 1) / n * sigma1_sq_hats[j]
  }
  
  # Funções para calcular Viés, EQM e Erros Padrão (SE)
  
  # Calcula Viés e SE do Viés para um vetor de estimativas (theta_hats)
  calcular_vies_se <- function(theta_hats) {
    # G_j = theta_hat - theta (para Bias)
    G_vies <- theta_hats - sigma2_real 
    Vies_hat <- mean(G_vies)
    
    # SE do Vies (SE(G_bar) = sd(G_j) / sqrt(M)) [5]
    SE_Vies <- sd(G_vies) / sqrt(M)
    return(c(Vies = Vies_hat, SE_Vies = SE_Vies))
  }
  
  # Calcula EQM e SE do EQM para um vetor de estimativas (theta_hats)
  calcular_eqm_se <- function(theta_hats) {
    # G_j = (theta_hat - theta)^2 (para MSE/EQM)
    G_eqm <- (theta_hats - sigma2_real)^2 
    EQM_hat <- mean(G_eqm)
    
    # SE do EQM (SE(G_bar) = sd(G_j) / sqrt(M)) [5]
    SE_EQM <- sd(G_eqm) / sqrt(M)
    return(c(EQM = EQM_hat, SE_EQM = SE_EQM))
  }
  
  # Resultados para Sigma1^2 (Não Viciado)
  vies1 <- calcular_vies_se(sigma1_sq_hats)
  eqm1  <- calcular_eqm_se(sigma1_sq_hats)
  
  # Resultados para Sigma2^2 (Viciado)
  vies2 <- calcular_vies_se(sigma2_sq_hats)
  eqm2  <- calcular_eqm_se(eqm2_hats)
  
  # Organizando em DataFrame
  res1 <- data.frame(
    Distribuicao = dist_name, n = n, Estimador = "Sigma1^2 (Não Viciado)",
    vies = vies1["Vies"], se_vies = vies1["SE_Vies"],
    eqm = eqm1["EQM"], se_eqm = eqm1["SE_EQM"]
  )
  res2 <- data.frame(
    Distribuicao = dist_name, n = n, Estimador = "Sigma2^2 (Viciado)",
    vies = vies2["Vies"], se_vies = vies2["SE_Vies"],
    eqm = eqm2["EQM"], se_eqm = eqm2["SE_EQM"]
  )
  
  return(rbind(res1, res2))
}

# 2. Estudo para Distribuição Normal N(10, 4), sigma^2 = 4
mu_norm <- 10
sigma2_norm <- 4
r_norm <- function(n) { rnorm(n, mean = mu_norm, sd = sqrt(sigma2_norm)) }

for (n in n_valores) {
  resultados[[length(resultados) + 1]] <- simulacao_variancia_mc(
    n = n, M = M, dist_name = "Normal (N(10, 4))", mu = mu_norm, 
    sigma2 = sigma2_norm, r_function = r_norm
  )
}

# 3. Estudo para Distribuição Exponencial Exp(0.5), sigma^2 = 1/0.5^2 = 4
lambda_exp <- 0.5
mu_exp <- 1/lambda_exp
sigma2_exp <- 1/lambda_exp^2 # Variância = 1/lambda^2 = 4
r_exp <- function(n) { rexp(n, rate = lambda_exp) }

for (n in n_valores) {
  resultados[[length(resultados) + 1]] <- simulacao_variancia_mc(
    n = n, M = M, dist_name = "Exponencial (Exp(0.5))", mu = mu_exp, 
    sigma2 = sigma2_exp, r_function = r_exp
  )
}

# 4. Apresentação dos resultados
tabela_final <- do.call(rbind, resultados)

# Arredondando para melhor leitura
tabela_final[, 4:7] <- round(tabela_final[, 4:7], 6)

print("Resultados do Estudo Monte Carlo (M=1000 réplicas): Viés e EQM da Variância")
print(tabela_final)