# Implementaremos os três métodos em R, replicando cada simulação M=200 vezes para avaliar a variabilidade empírica (erro padrão) de cada estimador.
# Usaremos K=10.000 réplicas em cada simulação para obter boas estimativas.
# Definição dos parâmetros
K <- 10000    # Número de amostras Monte Carlo por replicação
M <- 200      # Número de replicações para estimar a variabilidade
mu_Y <- 1     # Média de Y
sigma2_Y <- 1 # Variância de Y (sigma=1)
sigma2_X_dado_Y <- 4 # Variância condicional de X|Y (sigma=2)
theta_exact <- 0.5 # Valor analítico exato de P(X > 1)

# Vetores para armazenar as M estimativas de cada método
estimates_mc <- numeric(M)
estimates_ce <- numeric(M)
estimates_ce_cv <- numeric(M)

# Loop principal para as 200 replicações
set.seed(42) 

for (m in 1:M) {
  
  # Geração das amostras de Y (N(1, 1))
  Y <- rnorm(K, mean = mu_Y, sd = sqrt(sigma2_Y))
  
  # -----------------------------------------------------------------
  # Método a) Monte Carlo Simples (MC)
  # -----------------------------------------------------------------
  
  # Geração de X a partir de Y e do erro (erro ~ N(0, 4))
  epsilon <- rnorm(K, mean = 0, sd = sqrt(sigma2_X_dado_Y))
  X <- Y + epsilon
  
  # Estimativa: Média de I(X > 1)
  estimates_mc[m] <- mean(X > 1)
  
  # -----------------------------------------------------------------
  # Método b) Esperança Condicional (CE)
  # -----------------------------------------------------------------
  
  # Cálculo da função g(Y) = Phi((Y - 1) / 2)
  g_Y <- pnorm((Y - 1) / 2)
  
  # Estimativa: Média de g(Y)
  estimates_ce[m] <- mean(g_Y)
  
  # -----------------------------------------------------------------
  # Método c) Esperança Condicional com Variável de Controle Y (CE + CV)
  # -----------------------------------------------------------------
  
  # 1. Cálculo do coeficiente ótimo estimado (c_hat)
  # Cov(g(Y), Y) / Var(Y)
  
  # Covariância amostral (usamos cov/var para estimar c a partir dos K pares)
  c_hat <- cov(g_Y, Y) / var(Y)
  
  # 2. Cálculo da Estimativa CV: theta_ce - c_hat * (Y_bar - E[Y])
  Y_bar <- mean(Y)
  theta_ce_hat <- estimates_ce[m] # Já calculada acima
  
  estimates_ce_cv[m] <- theta_ce_hat - c_hat * (Y_bar - mu_Y)
}

# -----------------------------------------------------------------
# Análise de Variabilidade e Comparação
# -----------------------------------------------------------------

# Variabilidade (Erro Padrão Empírico das M=200 replicações)
se_mc <- sd(estimates_mc)
se_ce <- sd(estimates_ce)
se_ce_cv <- sd(estimates_ce_cv)

# Média das Estimativas (Viés empírico)
mean_mc <- mean(estimates_mc)
mean_ce <- mean(estimates_ce)
mean_ce_cv <- mean(estimates_ce_cv)

# Criação de um Data Frame para os resultados
results <- data.frame(
  Metodo = c("MC Simples (a)", "Esperança Condicional (b)", "CE + Variável de Controle (c)"),
  Media_Estimada = c(mean_mc, mean_ce, mean_ce_cv),
  Erro_Padrao_Estimado = c(se_mc, se_ce, se_ce_cv)
)

print(paste("Valor Exato de Theta (P{X > 1}):", theta_exact))
print(results)

# Cálculo da Redução Percentual na Variância
# Comparação com o MC Simples (a)
reduction_ce <- (se_mc^2 - se_ce^2) / se_mc^2 * 100
reduction_ce_cv <- (se_mc^2 - se_ce_cv^2) / se_mc^2 * 100

cat("\n--- Redução de Variância (Comparada ao MC Simples) ---\n")
cat("Redução de Variância (b) vs (a):", round(reduction_ce, 2), "%\n")
cat("Redução de Variância (c) vs (a):", round(reduction_ce_cv, 2), "%\n")

# -----------------------------------------------------------------
# Visualização da Variabilidade (d)
# -----------------------------------------------------------------

boxplot(estimates_mc, estimates_ce, estimates_ce_cv, 
        names = c("MC Simples", "Condicional", "Cond. + CV"),
        main = "Variabilidade dos Estimadores (M=200)",
        ylab = "Estimativas de Theta",
        col = c("lightblue", "lightgreen", "salmon"))
abline(h = theta_exact, col = "red", lty = 2, lwd = 2)