# 1. Função manual para a FDP da distribuição Normal (dnorm_manual)
# Esta função calcula phi(x | mu, sigma^2)
dnorm_manual <- function(x, mu, sigma) {
  if (sigma <= 0) stop("O desvio padrão (sigma) deve ser positivo.")
  
  exponent <- -((x - mu)^2 / (2 * sigma^2))
  coefficient <- 1 / (sigma * sqrt(2 * pi))
  return(coefficient * exp(exponent))
}

# 2. Função manual para a FDA da distribuição Normal (pnorm_manual)
# Esta função calcula Phi(q | mu, sigma^2) usando integrate
pnorm_manual <- function(q, mu, sigma) {
  if (sigma <= 0) stop("O desvio padrão (sigma) deve ser positivo.")
  
  # Casos extremos para q, onde a FDA é 0 ou 1
  if (q == -Inf) return(0)
  if (q == Inf) return(1)
  
  # Usar integrate para calcular a área sob a curva da FDP
  # desde -Inf até q.
  # A função integrate retorna uma lista, e o valor da integral está em 'value'.
  result <- tryCatch({
    integrate(f = function(x) dnorm_manual(x, mu, sigma),
              lower = -Inf, upper = q)
  }, error = function(e) {
    stop(paste("Erro durante a integração em pnorm_manual para q =", q, ":", e$message))
  })
  
  return(result$value)
}

# 3. Função reimplementada para gerar amostras da Normal Truncada
# Esta função calcula x = F_trunc_inv(u) usando pnorm_manual e qnorm
rnormtrunc_reimplemented <- function(n, mu, sigma, a, b) {
  # Validação de entradas
  if (sigma <= 0) stop("O desvio padrão (sigma) deve ser positivo.")
  if (a >= b) stop("O limite inferior 'a' deve ser menor que o limite superior 'b'.")
  
  # Gerar n números aleatórios uniformes entre 0 e 1
  us <- runif(n)
  
  # Calcular os valores da FDA não truncada nos pontos de truncamento 'a' e 'b'
  # Estes valores são constantes para todas as amostras e calculados uma vez.
  # Usamos max(0, min(1, ...)) para garantir que os valores estejam dentro de
  # para evitar erros com qnorm devido a pequenas imprecisões de integração.
  Phi_a <- max(0, min(1, pnorm_manual(a, mu, sigma)))
  Phi_b <- max(0, min(1, pnorm_manual(b, mu, sigma)))
  
  # Calcular P_u conforme a derivação matemática
  # Esta é a probabilidade acumulada ajustada que será usada em qnorm
  P_u <- Phi_a + us * (Phi_b - Phi_a)
  
  # Gerar as amostras da distribuição Normal truncada usando qnorm
  # qnorm(p, mean, sd) retorna o quantil para a probabilidade p
  amostra <- qnorm(P_u, mean = mu, sd = sigma)
  
  return(amostra)
}

# Função dnormtrunc original (do material fornecido) para comparação gráfica
# Esta função usa as funções dnorm e pnorm do R, que são permitidas para comparação
dnormtrunc <- function(x, mu, sigma, a, b){
  d <- dnorm(x, mu, sigma) / (pnorm(b, mu, sigma) - pnorm(a, mu, sigma))
  return(d)
}

### Exemplo de Uso e Ilustração

set.seed(123456789) # Para reprodutibilidade

# Parâmetros da distribuição Normal truncada
n_samples <- 10000
mu_val <- 10
sigma_val <- 2
a_val <- 6
b_val <- 13

# Gerar amostras usando a função reimplementada
amostra_reimplementada <- rnormtrunc_reimplemented(n_samples, mu_val, sigma_val, a_val, b_val)

# Configurações para o gráfico
par(mar = c(4, 4, 1, 1))

# Criar um histograma das amostras geradas
hist(amostra_reimplementada, xlab = "x", ylab = "Densidade", prob = TRUE, main = "",
     xlim = c(a_val, b_val), col = "lightblue", border = "darkblue",
     ylim = c(0, max(dnormtrunc(seq(a_val, b_val, length.out = 100), mu_val, sigma_val, a_val, b_val)) * 1.1))

# Adicionar a curva de densidade teórica para comparação
grid_x <- seq(a_val, b_val, length.out = 500)
lines(grid_x, dnormtrunc(grid_x, mu_val, sigma_val, a_val, b_val), col = "red", lwd = 2)

legend("topleft", legend = c("Amostras (Reimplementado)", "Densidade Teórica"),
       col = c("darkblue", "red"), lty = c(NA, 1), pch = c(22, NA), lwd = c(NA, 2),
       fill = c("lightblue", NA), border = c("darkblue", NA), bty = "n")