# Para fins de simulação de Monte Carlo, a integral θ pode ser vista como a esperança de g(U), onde g(U)=e U e U∼Uniforme(0,1).
# Usamos k=100.000 réplicas para garantir uma boa aproximação.

set.seed(123) 

# 1. Parâmetros da Simulação
k <- 100000 # Número de réplicas Monte Carlo
g <- function(x) {
  # Função a ser integrada: e^x
  return(exp(x))
}

# -----------------------------------------------
# 2. Método de Monte Carlo Simples (MC)
# -----------------------------------------------

# Geração de amostras U_i ~ Uniforme(0, 1) [2, 4]
U <- runif(k) 

# Variáveis Monte Carlo (X_i = g(U_i))
X_mc <- g(U)

# Estimativa de theta: média amostral [5, 6]
theta_mc_hat <- mean(X_mc)

# Estimativa da Variância de Var[g(U)]
# Usamos var(X_mc) que calcula a variância amostral com denominador k-1
var_g_u_mc <- var(X_mc) 

# Estimativa da Variância do Estimador MC: Var[theta_hat] = Var[g(U)] / k [7]
var_theta_mc_hat <- var_g_u_mc / k

# -----------------------------------------------
# 3. Método da Variável Antitética (AV)
# -----------------------------------------------

# A variável antitética é U_ant = 1 - U, que também é Uniforme(0, 1) [Baseado em conhecimento externo]
U_antithethic <- 1 - U

# Gerando as duas estimativas correlacionadas
X_ant_1 <- g(U)          # e^U
X_ant_2 <- g(U_antithethic) # e^(1-U)

# Variável antitética Y_i = (X_ant_1 + X_ant_2) / 2
Y_ant <- (X_ant_1 + X_ant_2) / 2

# Estimativa de theta: média amostral de Y_ant [6]
theta_ant_hat <- mean(Y_ant)

# Estimativa da Variância de Var[Y_ant]
var_Y_ant <- var(Y_ant)

# Estimativa da Variância do Estimador AV: Var[theta_ant_hat] = Var[Y_ant] / k [7]
var_theta_ant_hat <- var_Y_ant / k

# -----------------------------------------------
# 4. Cálculo da Redução Percentual na Variância
# -----------------------------------------------

# Redução percentual: (Var_MC - Var_ANT) / Var_MC * 100
reduction_pct <- ((var_theta_mc_hat - var_theta_ant_hat) / var_theta_mc_hat) * 100

# -----------------------------------------------
# 5. Apresentação dos Resultados
# -----------------------------------------------

theta_exact <- exp(1) - 1

cat("=========================================================\n")
cat("Análise da Estimativa de θ = ∫[1] exp(x) dx \n")
cat("Valor Analítico Exato (e - 1):", round(theta_exact, 7), "\n")
cat("Número de Réplicas (k):", k, "\n")
cat("=========================================================\n")

cat("\n--- Monte Carlo Simples (MC) ---\n")
cat("Estimativa (θ_hat_MC):", round(theta_mc_hat, 7), "\n")
cat("Variância Empírica Estimada:", format(var_theta_mc_hat, scientific = TRUE, digits = 5), "\n")
cat("Erro Padrão Empírico (SE):", format(sqrt(var_theta_mc_hat), scientific = TRUE, digits = 5), "\n")

cat("\n--- Variável Antitética (AV) ---\n")
cat("Estimativa (θ_hat_ANT):", round(theta_ant_hat, 7), "\n")
cat("Variância Empírica Estimada:", format(var_theta_ant_hat, scientific = TRUE, digits = 5), "\n")
cat("Erro Padrão Empírico (SE):", format(sqrt(var_theta_ant_hat), scientific = TRUE, digits = 5), "\n")

cat("\n--- Redução de Variância ---\n")
cat("Redução Percentual na Variância:", round(reduction_pct, 4), "%\n")
cat("=========================================================\n")