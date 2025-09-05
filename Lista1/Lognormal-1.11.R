# Função para gerar variáveis aleatórias Normais padrão (N(0,1))
rnormal_padrao <- function(n){
  # Calcula o número de iterações necessárias para obter n amostras.
  # Cada iteração bem-sucedida gera 2 valores normais padrão.
  itera <- ceiling(n/2) 
  amostra <- numeric(itera * 2) # Pré-aloca espaço para a amostra
  
  for(i in 1:itera){
    s <- 2 # Inicializa s com um valor > 1 para entrar no loop while
    while(s > 1 || s == 0){ # A condição s == 0 é adicionada para robustez, embora rara
      u1 <- runif(1)
      u2 <- runif(1)
      v1 <- 2*u1-1 # Transforma U(0,1) para U(-1,1)
      v2 <- 2*u2-1 # Transforma U(0,1) para U(-1,1)
      s <- v1^2 + v2^2 # Calcula s = R^2
    }
    
    # Aplica as transformações de Box-Muller para obter X e Y ~ N(0,1)
    scale_factor <- sqrt(-2*log(s)/s)
    x <- scale_factor * v1
    y <- scale_factor * v2
    
    # Armazena os valores gerados
    amostra[i*2-1] <- x
    amostra[i*2] <- y
  }
  return(amostra[1:n]) # Retorna exatamente n amostras
}

# Função para gerar variáveis aleatórias Lognormais usando o método da transformação
rlnorm_transformacao <- function(n, meanlog = 1, sdlog = 1) {
  # Verifica se sdlog é positivo para uma distribuição não degenerada
  if (sdlog <= 0) {
    warning("sdlog deve ser positivo para uma distribuição Lognormal não degenerada. Retornando valores constantes.")
    return(rep(exp(meanlog), n))
  }
  
  # 1. Gerar variáveis Normais padrão (Z_std ~ N(0,1))
  z_std <- rnormal_padrao(n)
  
  # 2. Transformar para Normal com os parâmetros desejados (Z ~ N(meanlog, sdlog^2))
  # Se Z ~ N(0,1), então mu + sigma*Z ~ N(mu, sigma^2)
  z_general <- meanlog + sdlog * z_std
  
  # 3. Exponenciar para obter Lognormal (X = e^Z)
  x_lognormal <- exp(z_general)
  
  return(x_lognormal)
}

### Método da Composição para Geração Lognormal

### Comparação com Histogramas

# Parâmetros da distribuição Lognormal
meanlog_param <- 1
sdlog_param <- 1
n_amostras <- 1000

# 1. Gerar amostras usando a função implementada (método da transformação)
set.seed(123) # Para reprodutibilidade
amostra_transformacao <- rlnorm_transformacao(n_amostras, meanlog = meanlog_param, sdlog = sdlog_param)

# 2. Gerar amostras usando a função dlnorm do R para comparação
amostra_r_dlnorm <- rlnorm(n_amostras, meanlog = meanlog_param, sdlog = sdlog_param)

# Definir o grid para a curva de densidade teórica
min_val <- min(amostra_transformacao, amostra_r_dlnorm)
max_val <- max(amostra_transformacao, amostra_r_dlnorm)
grid_x <- seq(min_val, max_val, length.out = 500)
densidade_teorica <- dlnorm(grid_x, meanlog = meanlog_param, sdlog = sdlog_param)

# Configurar o layout para dois histogramas lado a lado
par(mfrow = c(1, 2), mar = c(5, 4, 4, 2) + 0.1)

# Histograma da amostra gerada pelo método da transformação
hist(amostra_transformacao, freq = FALSE, breaks = 30,
     main = "Lognormal (Transformação)",
     xlab = "Valor", ylab = "Densidade",
     col = "lightblue", border = "darkblue",
     xlim = c(min_val, max_val), ylim = c(0, max(densidade_teorica) * 1.1))
lines(grid_x, densidade_teorica, col = "red", lwd = 2)
legend("topleft", legend = c("Densidade Teórica"), col = "red", lty = 1, lwd = 2, bty = "n")

# Histograma da amostra gerada pela função rlnorm do R
hist(amostra_r_dlnorm, freq = FALSE, breaks = 30,
     main = "Lognormal (rlnorm do R)",
     xlab = "Valor", ylab = "Densidade",
     col = "lightgreen", border = "darkgreen",
     xlim = c(min_val, max_val), ylim = c(0, max(densidade_teorica) * 1.1))
lines(grid_x, densidade_teorica, col = "red", lwd = 2)
legend("topleft", legend = c("Densidade Teórica"), col = "red", lty = 1, lwd = 2, bty = "n")