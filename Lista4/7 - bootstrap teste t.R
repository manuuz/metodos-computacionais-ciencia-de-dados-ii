# --- 1. Funcao do Teste Bootstrap de Media (Baseado em 3.6.2) ---
teste_boot_1mean <- function(z, mu0, B = 1000) {
  n <- length(z)
  
  # 1. Estatística t observada (t_obs)
  # Usamos a estatística t do teste padrão para a amostra original
  teste_t_obs <- t.test(z, mu = mu0)
  t_obs <- teste_t_obs$statistic
  
  # 2. Gerar amostras Bootstrap centradas sob H0 (3.6.2)
  # A distribuicao empirica e transladada para ter media mu0.
  z_til <- z - mean(z) + mu0
  
  t_ast <- numeric(B)
  
  for(b in 1:B){
    # Amostra bootstrap com reposicao da amostra centrada
    z_ast <- sample(z_til, size = n, replace = TRUE)
    
    # Calculamos a estatistica t* na amostra centrada para encontrar a distribuicao sob H0
    # t* = (media_z_ast - mu0) / (se_z_ast)
    
    # O teste t.test em R calcula o denominador da estatística como sd(z_ast) / sqrt(n)
    # Note que mean(z_ast) deveria ser proximo de mu0, mas usamos a expressao de t.test
    
    # Usamos o resultado do t.test, que e a estatistica t-studentizada
    t_ast[b] <- t.test(z_ast, mu = mu0)$statistic
  }
  
  # 3. Calculo do p-valor bootstrap
  # O p-valor e a proporcao de estatisticas t* que sao mais extremas que t_obs (3.6.2)
  p_valor <- ( sum(t_ast >= abs(t_obs)) + sum(t_ast <= -abs(t_obs)) ) / B
  
  return(p_valor)
}

# --- 2. Funcao Principal do Estudo Monte Carlo para Poder ---
estudo_poder_mc <- function(mu_vec, n, M, B, alpha) {
  
  num_mu <- length(mu_vec)
  poder_t <- numeric(num_mu)
  poder_boot <- numeric(num_mu)
  
  # O teste e H0: mu = 2 (mu0 = 2)
  mu0 <- 2
  
  for (i in 1:num_mu) {
    mu_real <- mu_vec[i]
    # Em uma distribuicao Chi-quadrado, E[X] = nu. Entao, nu = mu_real.
    nu_real <- mu_real 
    
    rejeicoes_t <- 0
    rejeicoes_boot <- 0
    
    # Loop Monte Carlo (M replicacoes)
    for (m in 1:M) {
      
      # 1. Geracao da Amostra (violando a normalidade)
      # Gerar dados da Chi-quadrado com E[X] = mu_real = nu_real
      dados_mc <- rchisq(n, df = nu_real)
      
      # --- A. Teste t padrao ---
      # O teste t e o teste de comparacao, que tem sua hipotese de normalidade violada.
      teste_t <- t.test(dados_mc, mu = mu0, alternative = "two.sided")
      if (teste_t$p.value <= alpha) {
        rejeicoes_t <- rejeicoes_t + 1
      }
      
      # --- B. Teste Bootstrap ---
      p_valor_boot <- teste_boot_1mean(dados_mc, mu0 = mu0, B = B)
      if (p_valor_boot <= alpha) {
        rejeicoes_boot <- rejeicoes_boot + 1
      }
    }
    
    # Estimativa da Funcao Poder (2.3)
    poder_t[i] <- rejeicoes_t / M
    poder_boot[i] <- rejeicoes_boot / M
  }
  
  return(data.frame(
    mu_real = mu_vec,
    Poder_T = poder_t,
    Poder_Bootstrap = poder_boot
  ))
}

# --- 3. Execucao do Estudo Monte Carlo ---

# Valores de mu (nu) para avaliar a funcao poder, incluindo o valor nulo mu=2
mu_valores <- c(1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0)
N_amostral <- 30
M_mc <- 5000
B_boot <- 1000
Alpha_nivel <- 0.05

set.seed(42)
resultados_poder <- estudo_poder_mc(
  mu_vec = mu_valores,
  n = N_amostral,
  M = M_mc,
  B = B_boot,
  alpha = Alpha_nivel
)

# Imprimir resultados
print("Resultados da Funcao Poder (M=5000, n=30, Populacao: Chi-quadrado)")
print(round(resultados_poder, 4))

# --- 4. Plotagem da Funcao Poder ---
library(ggplot2)

df_plot <- data.frame(
  mu_real = resultados_poder$mu_real,
  Poder = c(resultados_poder$Poder_T, resultados_poder$Poder_Bootstrap),
  Metodo = factor(rep(c("Teste T (Violado)", "Bootstrap"), each = length(mu_valores)))
)

p_power <- ggplot(df_plot, aes(x = mu_real, y = Poder, color = Metodo)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_vline(xintercept = 2, linetype = "dashed", color = "black", alpha = 0.7) +
  geom_hline(yintercept = Alpha_nivel, linetype = "dotted", color = "gray50") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Função Poder Estimada: Teste T vs. Bootstrap",
    subtitle = "População: Qui-quadrado (Assimetria viola a suposição de Normalidade)",
    x = expression(paste("Média Populacional Real (", mu, ")")),
    y = "Poder (P(Rejeitar H0))"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p_power)