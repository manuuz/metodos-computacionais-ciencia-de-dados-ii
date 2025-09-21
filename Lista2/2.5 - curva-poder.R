# ==============================================================================
# EXERCÍCIO 2.5 (Adaptado do Exemplo 3.7/2.7) - CURVA DE PODER BILATERAL
# ==============================================================================

# Definindo parâmetros do estudo
M <- 5000           # Número de réplicas Monte Carlo para estimar o poder
alpha <- 0.05      # Nível de significância
mu0 <- 500         # Média sob H0
sigma <- 100       # Desvio padrão populacional (assumido do Exemplo 2.7) [1]
n_valores <- c(10, 20, 30, 40, 50) # Tamanhos amostrais solicitados

# Definindo o intervalo de médias sob a hipótese alternativa (mu_real)
# Usamos um intervalo simétrico em torno de mu0 para plotar a curva de poder.
mu_alternativa <- seq(400, 600, by = 5)

# Matriz para armazenar as estimativas de poder (linhas = mu, colunas = n)
power_matrix <- matrix(NA, nrow = length(mu_alternativa), ncol = length(n_valores))
colnames(power_matrix) <- paste0("n=", n_valores)

set.seed(42) # Para reprodutibilidade

# Loop 1: Iterar sobre os tamanhos amostrais (n)
for (j in 1:length(n_valores)) {
  n <- n_valores[j]
  
  # Loop 2: Iterar sobre os valores reais da média (mu_real)
  for (i in 1:length(mu_alternativa)) {
    mu_real <- mu_alternativa[i]
    
    # Simulação MC para estimar o poder pi(mu_real) [3]
    pvalues <- replicate(M, expr = {
      # 1. Gerar amostra sob mu_real (assumindo Normal, como no Exemplo 2.7) [1]
      x <- rnorm(n, mean = mu_real, sd = sigma) 
      
      # 2. Executar o teste t BILATERAL (alternative = "two.sided") [4]
      teste_t <- t.test(x, alternative = "two.sided", mu = mu0)
      
      # Retorna o p-valor
      teste_t$p.value
    })
    
    # 3. Estimar o poder (Proporção de rejeições H0: p-value < alpha) [1]
    power_matrix[i, j] <- mean(pvalues < alpha)
  }
}

# ==============================================================================
# PLOTAGEM DA CURVA DE PODER EMPÍRICA
# ==============================================================================

# Definindo cores e legendas
cores <- c("blue", "red", "green4", "purple", "orange")
nomes_curvas <- paste("n =", n_valores)

par(mar = c(5, 4, 4, 2) + 0.1)

# Plot inicial da primeira curva (n=10)
plot(mu_alternativa, power_matrix[, 1], 
     type = "l", 
     col = cores[5], 
     lwd = 2,
     ylim = c(0, 1),
     xlab = expression(paste("Média Populacional Real (", mu, ")")), 
     ylab = "Poder (Probabilidade de Rejeição)",
     main = "Curva de Poder Empírica para o Teste t Bilateral",
     las = 1)

# Adicionar as curvas restantes (n=20 a n=50)
for (j in 2:length(n_valores)) {
  lines(mu_alternativa, power_matrix[, j], col = cores[j], lwd = 2)
}

# Linhas de referência: Nível de significância (h=0.05) e Média Nula (v=500)
abline(h = alpha, lty = 2, col = "gray50")
abline(v = mu0, lty = 3, col = "gray50")

# Adicionar legenda
legend("topright", 
       legend = nomes_curvas, 
       col = cores, 
       lwd = 2,
       title = "Tamanho Amostral (n)",
       bty = "n")

print("Poder Estimado (Matriz de Resultados):")
print(round(power_matrix[c(1, 11, 21, 31, 41), ], 4))