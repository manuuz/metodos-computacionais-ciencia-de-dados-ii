# ==============================================================================
# FUNÇÃO 1: TESTE DE HIPÓTESES MONTE CARLO (THMC)
# ==============================================================================
# Testa H0: mu = mu0 vs H1: mu != mu0 (bilateral)
# Baseado na estatística Z, mas p-valor calculado via simulação.

thmc_ztest <- function(x, mu0, sigma, B = 10000) {
  # x: amostra observada
  # mu0: valor sob a hipótese nula
  # sigma: desvio padrão populacional (conhecido)
  # B: número de réplicas Monte Carlo para estimar o p-valor
  
  n <- length(x)
  se_mu <- sigma / sqrt(n)
  
  # 1. Calcular o valor observado da estatística Z
  z_obs <- (mean(x) - mu0) / se_mu
  
  # 2. Simular a estatística Z sob H0 (mu = mu0) B vezes
  # Geramos B amostras e calculamos suas médias (x_bar*)
  # Assumimos X ~ N(mu0, sigma) sob H0
  
  # Gerar B médias amostrais sob H0: x_bar* ~ N(mu0, sigma/sqrt(n))
  x_bar_star <- rnorm(B, mean = mu0, sd = se_mu)
  
  # Calcular a estatística Z* para cada média simulada
  z_star <- (x_bar_star - mu0) / se_mu
  
  # 3. Calcular o p-valor MC (proporção de |Z*| >= |Z_obs|) [3]
  p_valor_mc <- sum(abs(z_star) >= abs(z_obs)) / B
  
  return(list(
    z_obs = z_obs,
    p_valor_mc = p_valor_mc
  ))
}

# ==============================================================================
# FUNÇÃO 2: TESTE NORMAL PADRÃO (TESTE Z ANALÍTICO)
# Implementado para comparação direta.
# ==============================================================================
z_test_analitico <- function(x, mu0, sigma) {
  n <- length(x)
  se_mu <- sigma / sqrt(n)
  
  z_obs <- (mean(x) - mu0) / se_mu
  
  # p-valor bilateral P(|Z| >= |z_obs|)
  p_valor_norm <- 2 * pnorm(-abs(z_obs))
  
  return(list(
    z_obs = z_obs,
    p_valor_norm = p_valor_norm
  ))
}

# ==============================================================================
# ESTUDO MC PARA COMPARAR O PODER DE TESTE
# ==============================================================================

# Definindo parâmetros do estudo
n <- 30            # Tamanho da amostra
sigma <- 5         # Desvio padrão populacional (conhecido)
mu0 <- 10          # Média sob H0
mu1 <- 11.5        # Média real (sob H1)
alfa <- 0.05       # Nível de significância [4]
M <- 10000         # Número de réplicas para o estudo de poder [1]
B_mc <- 5000       # Número de réplicas internas para o THMC

# Vetores para armazenar as decisões de rejeição
rejeicoes_mc <- numeric(M)
rejeicoes_norm <- numeric(M)

set.seed(42) # Para reprodutibilidade

# Início do loop Monte Carlo para estimar o poder
for (i in 1:M) {
  
  # Passo 1: Gerar a amostra sob a hipótese alternativa (mu = mu1) [1]
  x_amostra <- rnorm(n, mean = mu1, sd = sigma)
  
  # --- THMC (Função 1) ---
  resultado_mc <- thmc_ztest(x_amostra, mu0 = mu0, sigma = sigma, B = B_mc)
  
  # Passo 3: Decisão de rejeição para o THMC
  if (resultado_mc$p_valor_mc < alfa) {
    rejeicoes_mc[i] <- 1 # Rejeita H0
  }
  
  # --- Teste Normal Analítico (Função 2) ---
  resultado_norm <- z_test_analitico(x_amostra, mu0 = mu0, sigma = sigma)
  
  # Passo 3: Decisão de rejeição para o Teste Normal
  if (resultado_norm$p_valor_norm < alfa) {
    rejeicoes_norm[i] <- 1 # Rejeita H0
  }
}

# 4. Estimar o poder de teste [2]
poder_mc <- mean(rejeicoes_mc)
poder_norm <- mean(rejeicoes_norm)

# 5. Estimar o erro padrão Monte Carlo para as estimativas de poder
# Fórmula do erro padrão para uma proporção: sqrt(p*(1-p)/M) [5]
se_poder_mc <- sqrt(poder_mc * (1 - poder_mc) / M)
se_poder_norm <- sqrt(poder_norm * (1 - poder_norm) / M)

# ==============================================================================
# RESULTADOS
# ==============================================================================

print("=== Comparação do Poder de Teste (mu real = 11.5) ===")

tabela_comparacao <- data.frame(
  Teste = c("THMC (Monte Carlo)", "Normal Analítico (Z-Test)"),
  Poder_Estimado = c(poder_mc, poder_norm),
  Erro_Padrao_MC = c(se_poder_mc, se_poder_norm)
)

print(tabela_comparacao)

# Análise de Convergência (Nível de Signif.: Testando mu=10 quando mu real=10)
# Para garantir que o THMC tem o nível nominal correto (0.05)
M_alfa <- 10000
rejeicoes_alfa <- numeric(M_alfa)
for (i in 1:M_alfa) {
  x_amostra <- rnorm(n, mean = mu0, sd = sigma)
  p_mc <- thmc_ztest(x_amostra, mu0 = mu0, sigma = sigma, B = B_mc)$p_valor_mc
  if (p_mc < alfa) {
    rejeicoes_alfa[i] <- 1
  }
}
alfa_empirico_mc <- mean(rejeicoes_alfa)
se_alfa_mc <- sqrt(alfa_empirico_mc * (1 - alfa_empirico_mc) / M_alfa)

print(sprintf("\n--- Nível de Significância Empírico (mu real = 10) ---"))
print(sprintf("Taxa empírica de Erro Tipo I (THMC): %.4f (SE: %.5f)", 
              alfa_empirico_mc, se_alfa_mc))