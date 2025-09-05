# Função para gerar amostras de uma distribuição categórica (rcategorica2 da fonte)
rcategorica2 <- function(n, x, p){
  # Validação de entradas
  if (length(x) != length(p)) stop("Os vetores de valores e probabilidades devem ter o mesmo comprimento.")
  if (sum(p) != 1) {
    warning("As probabilidades não somam 1. Elas serão normalizadas.")
    p <- p / sum(p)
  }
  if (any(p < 0)) stop("As probabilidades devem ser não negativas.")
  
  # Ordenar os valores e probabilidades em ordem decrescente de probabilidade
  ordem <- order(p, decreasing = TRUE)
  x_ordered <- x[ordem]
  p_ordered <- p[ordem]
  
  am <- numeric(n) # Vetor para armazenar as amostras
  acum <- cumsum(p_ordered) # Calcular as probabilidades acumuladas
  
  # Gerar n amostras
  for(i in 1:n){
    u <- runif(1) # Gerar um número aleatório uniforme
    # Encontrar a categoria correspondente
    for(j in 1:length(x_ordered)){
      if(u <= acum[j]){
        am[i] <- x_ordered[j]
        break # Sair do loop interno após encontrar a categoria
      }
    }
  }
  return(am)
}

set.seed(123) # Para reprodutibilidade

# Número de amostras para cada cenário de simulação
n_samples_sim <- 50000

# Função auxiliar para plotar comparações (adaptada da ilustração da fonte)
plot_comparison <- function(amostra, theoretical_p, title_str, x_vals) {
  # Assegurar que os níveis do fator correspondam aos valores possíveis para ordenação correta
  observed_freq <- prop.table(table(factor(amostra, levels = x_vals)))
  
  # Encontrar o máximo valor para o limite do eixo Y
  max_y <- max(max(observed_freq), max(theoretical_p)) * 1.1
  
  # Plotar o histograma das frequências observadas
  barplot(observed_freq, ylab = "Frequência Relativa", xlab = "Valor de X",
          main = title_str, col = rainbow(length(x_vals), s = 0.7, v = 0.9), ylim = c(0, max_y))
  
  # Adicionar as probabilidades teóricas como pontos para comparação
  points(x = seq_along(x_vals), y = theoretical_p, col = "red", pch = 16, cex = 1.5)
  
  # Adicionar legenda
  legend("topright", legend = c("Frequência Observada", "Probabilidade Teórica"),
         col = c(NA, "red"), pch = c(NA, 16), lty = c(1, NA), bty = "n",
         fill = c(rainbow(length(x_vals), s = 0.7, v = 0.9), NA), border = c("black", NA))
}

# Cenário 1: Distribuição Uniforme Discreta (quatro categorias)
cat_x_1 <- 1:4
cat_p_1 <- c(0.25, 0.25, 0.25, 0.25)
sim_amostra_1 <- rcategorica2(n_samples_sim, cat_x_1, cat_p_1)

# Cenário 2: Distribuição Categórica Assimétrica (valores numéricos distintos)
cat_x_2 <- c(10, 20, 30, 40, 50)
cat_p_2 <- c(0.05, 0.15, 0.30, 0.25, 0.25)
sim_amostra_2 <- rcategorica2(n_samples_sim, cat_x_2, cat_p_2)

# Cenário 3: Distribuição Categórica com muitos níveis e valores não numéricos
cat_x_3 <- LETTERS[1:7] # Usando letras como possíveis valores
cat_p_3 <- c(0.05, 0.1, 0.15, 0.2, 0.15, 0.2, 0.15)
sim_amostra_3 <- rcategorica2(n_samples_sim, cat_x_3, cat_p_3)

# Plotar os resultados dos três cenários
par(mfrow = c(3, 1), mar = c(4, 4, 3, 1))

plot_comparison(sim_amostra_1, cat_p_1, "Cenário 1: Distribuição Uniforme Discreta", cat_x_1)
plot_comparison(sim_amostra_2, cat_p_2, "Cenário 2: Distribuição Categórica Assimétrica (Numérica)", cat_x_2)
plot_comparison(sim_amostra_3, cat_p_3, "Cenário 3: Distribuição Categórica com Caracteres", cat_x_3)

# Restaurar as configurações gráficas
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)