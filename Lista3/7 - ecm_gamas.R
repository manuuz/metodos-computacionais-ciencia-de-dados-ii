# Esta implementação em R segue as formulações do algoritmo EM/ECM para a mistura de duas distribuições Gama, conforme detalhado nas fontes. O algoritmo é considerado ECM (Expectation-Conditional Maximization) porque a maximização dos parâmetros de forma requer a solução numérica de equações não lineares.
# 1. Funções Auxiliares
# Precisamos da função de densidade Gama, da função de probabilidade p (Passo E) e de uma função para resolver numericamente os parâmetros de forma, que envolve a função digama.
# Requer a biblioteca 'numDeriv' para calcular a função digamma (embora 'digamma' seja nativa no R) e 'stats' para otimização.
library(stats) 

# Função de densidade Gama usada no modelo (assumindo parametrização shape/rate)
dgamma_mix <- function(x, alpha, beta) {
  # f(xi|alpha_k, beta_k) = (x^(alpha_k - 1) * exp(-x * beta_k) * (beta_k^alpha_k)) / gamma(alpha_k) [1, 7]
  # Utilizamos a função nativa do R: dgamma(x, shape = alpha, rate = beta)
  return(dgamma(x, shape = alpha, rate = beta))
}

# Função para resolver numericamente a equação não linear para alpha (CM Step para alpha)
# Resolve Psi(alpha) = C_k, onde C_k é a média ponderada das log-observações mais log(beta_k)
# A solução deve ser positiva (alpha > 0)
solve_alpha <- function(C_k, alpha_init) {
  # Função cuja raiz deve ser encontrada: Psi(alpha) - C_k = 0
  H_k <- function(alpha) {
    if (alpha <= 0) return(NA)
    return(digamma(alpha) - C_k)
  }
  
  # Usamos uniroot para encontrar a raiz
  # O intervalo de busca precisa ser definido. Como Psi(alpha) cresce monotonicamente,
  # podemos encontrar limites razoáveis (por exemplo, usando a inversa aproximada da digamma)
  
  # Estimativa inicial rápida para definir o intervalo (opcional, mas recomendado para estabilidade)
  # Se C_k for grande, alpha é grande; se C_k for pequeno (negativo), alpha é pequeno.
  if (C_k > 20) { # Chute grande
    low = C_k - 5 
    high = C_k + 5 
  } else if (C_k > 0) { # Chute na faixa mais comum (1 a 10)
    low = 0.01
    high = 100
  } else { # Chute para alphas pequenos
    low = 0.001
    high = 5
  }
  
  # Se o intervalo não contiver uma raiz, tentamos um chute mais amplo
  if (sign(H_k(low)) == sign(H_k(high)) || is.na(H_k(low)) || is.na(H_k(high))) {
    low = 0.0001
    high = 500
  }
  
  # Garantindo que a solução numérica não falhe
  result <- try(uniroot(H_k, interval = c(low, high), tol = 1e-6)$root, silent = TRUE)
  
  if (inherits(result, "try-error")) {
    # Em caso de falha de uniroot (geralmente convergência ou intervalo incorreto), retorna a estimativa inicial
    return(alpha_init)
  }
  return(result)
}

### 2. Algoritmo ECM para Mistura de Duas Gamas

ECM_Gamma <- function(X, initial_params, max_iter = 500, tol = 1e-5) {
  
  n <- length(X)
  
  # Inicialização: theta = (delta, alpha1, alpha2, beta1, beta2)
  theta <- initial_params
  
  # Vetores para armazenar a convergência
  history <- data.frame(delta=numeric(max_iter), alpha1=numeric(max_iter), 
                        alpha2=numeric(max_iter), beta1=numeric(max_iter), 
                        beta2=numeric(max_iter))
  
  for (t in 1:max_iter) {
    
    delta_t <- theta[8]
    alpha1_t <- theta[9]
    alpha2_t <- theta[10]
    beta1_t <- theta[11]
    beta2_t <- theta[12]
    
    # ---------------------------------------------
    # PASSO E (Expectation Step)
    # ---------------------------------------------
    
    # Calculando as densidades condicionais: f(x_i | alpha_k, beta_k)
    f1 <- dgamma_mix(X, alpha1_t, beta1_t)
    f2 <- dgamma_mix(X, alpha2_t, beta2_t)
    
    # Calculando a probabilidade posterior (E[Z_i | x_i, theta^(t)] = p_i) [2, 13]
    numerator <- delta_t * f1
    denominator <- numerator + (1 - delta_t) * f2
    
    # Tratamento para evitar divisão por zero ou NA's
    p_i <- ifelse(denominator == 0, 0, numerator / denominator)
    
    # Estatísticas ponderadas
    sum_p <- sum(p_i)
    sum_1_minus_p <- n - sum_p # sum(1 - p_i)
    
    sum_x_p <- sum(X * p_i)
    sum_x_1_minus_p <- sum(X) - sum_x_p # sum(X * (1 - p_i))
    
    
    # ---------------------------------------------
    # PASSO M/CM (Maximization/Conditional Maximization Steps)
    # Sequência de maximizações para obter theta^(t+1)
    # ---------------------------------------------
    
    # 1. Atualização Analítica para delta [2, 13]
    # delta^(t+1) = Sum(p_i) / n
    delta_new <- sum_p / n
    
    # 2. Atualizações Analíticas para beta_k (dependem de alpha_k^(t)) [2, 6]
    
    # Para o Grupo 1
    if (sum_p > 0) {
      beta1_new <- alpha1_t * sum_p / sum_x_p
    } else {
      beta1_new <- beta1_t
    }
    
    # Para o Grupo 2
    if (sum_1_minus_p > 0) {
      beta2_new <- alpha2_t * sum_1_minus_p / sum_x_1_minus_p
    } else {
      beta2_new <- beta2_t
    }
    
    # 3. Atualizações Numéricas para alpha_k (Passos CM, usam beta_k^(t+1) fixo) [3, 6]
    
    # Para o Grupo 1: Resolver Psi(alpha1) = C1
    if (sum_p > 0) {
      C1 <- (sum(p_i * log(X)) / sum_p) + log(beta1_new)
      alpha1_new <- solve_alpha(C1, alpha1_t)
    } else {
      alpha1_new <- alpha1_t
    }
    
    # Para o Grupo 2: Resolver Psi(alpha2) = C2
    if (sum_1_minus_p > 0) {
      C2 <- (sum((1 - p_i) * log(X)) / sum_1_minus_p) + log(beta2_new)
      alpha2_new <- solve_alpha(C2, alpha2_t)
    } else {
      alpha2_new <- alpha2_t
    }
    
    # Novo vetor de parâmetros
    theta_new <- c(delta_new, alpha1_new, alpha2_new, beta1_new, beta2_new)
    
    # ---------------------------------------------
    # CRITÉRIO DE CONVERGÊNCIA
    # ---------------------------------------------
    
    diff_abs <- abs(theta_new - theta)
    denom <- abs(theta) + tol # Evita divisão por zero se o parâmetro for pequeno
    crit <- max(diff_abs / denom)
    
    theta <- theta_new
    history[t, ] <- theta
    
    if (crit < tol) {
      # Convergência alcançada
      history <- history[1:t, ]
      break
    }
  }
  
  if (t == max_iter && crit >= tol) {
    warning("O algoritmo ECM não convergiu no número máximo de iterações.")
  }
  
  # Renomeando as colunas de saída para maior clareza
  colnames(history) <- c("delta", "alpha1", "alpha2", "beta1", "beta2")
  
  return(list(
    estimates = theta,
    iterations = t,
    history = history
  ))
}

### 3. Geração de Dados e Avaliação

# Para avaliar o algoritmo, geramos dados a partir do modelo de mistura de duas Gamas.

# Função para gerar amostras da mistura de Gamas
rgamma_mixture <- function(n, delta, alpha1, alpha2, beta1, beta2) {
  # 1. Gerar os indicadores de grupo (Variável latente Z) [15]
  Z <- rbinom(n, size = 1, prob = delta)
  
  # 2. Gerar as amostras: Gama 1 se Z=1, Gama 2 se Z=0
  X <- numeric(n)
  
  # Grupo 1
  idx1 <- which(Z == 1)
  if (length(idx1) > 0) {
    X[idx1] <- rgamma(length(idx1), shape = alpha1, rate = beta1)
  }
  
  # Grupo 2
  idx2 <- which(Z == 0)
  if (length(idx2) > 0) {
    X[idx2] <- rgamma(length(idx2), shape = alpha2, rate = beta2)
  }
  
  return(X)
}

# ------------------------------------------------------------------------
# CENÁRIO 1: Mistura bem separada (Parâmetros de escala e forma diferentes)
# ------------------------------------------------------------------------

# Parâmetros Verdadeiros
TRUE_THETA1 <- c(delta = 0.6, alpha1 = 5, alpha2 = 2, beta1 = 1, beta2 = 0.2)
N <- 500
set.seed(42)
Data1 <- rgamma_mixture(N, TRUE_THETA1[8], TRUE_THETA1[9], TRUE_THETA1[10], TRUE_THETA1[11], TRUE_THETA1[12])

# Estimativas Iniciais (Chute longe, mas razoável)
# delta=0.5, alpha1=1, alpha2=1, beta1=0.5, beta2=0.5
INITIAL_THETA1 <- c(0.5, 1, 1, 0.5, 0.5)

cat("\n--- Cenário 1: Mistura bem separada (N=500) ---\n")
ECM_Result1 <- ECM_Gamma(Data1, INITIAL_THETA1)
print(paste("Iterações para convergência:", ECM_Result1$iterations))
df_estimates1 <- data.frame(
  Parameter = c("delta", "alpha1", "alpha2", "beta1", "beta2"),
  True = TRUE_THETA1,
  Estimate = ECM_Result1$estimates
)
print(round(df_estimates1, 5))

# ------------------------------------------------------------------------
# CENÁRIO 2: Mistura com tamanhos desiguais (Maior peso no Grupo 2)
# ------------------------------------------------------------------------

TRUE_THETA2 <- c(delta = 0.2, alpha1 = 10, alpha2 = 5, beta1 = 2, beta2 = 0.5)
N <- 1000
set.seed(100)
Data2 <- rgamma_mixture(N, TRUE_THETA2[8], TRUE_THETA2[9], TRUE_THETA2[10], TRUE_THETA2[11], TRUE_THETA2[12])
INITIAL_THETA2 <- c(0.5, 1, 1, 0.5, 0.5)

cat("\n--- Cenário 2: Mistura desigual (N=1000) ---\n")
ECM_Result2 <- ECM_Gamma(Data2, INITIAL_THETA2)
print(paste("Iterações para convergência:", ECM_Result2$iterations))
df_estimates2 <- data.frame(
  Parameter = c("delta", "alpha1", "alpha2", "beta1", "beta2"),
  True = TRUE_THETA2,
  Estimate = ECM_Result2$estimates
)
print(round(df_estimates2, 5))

# ------------------------------------------------------------------------
# CENÁRIO 3: Mistura mais difícil (densidades mais sobrepostas)
# ------------------------------------------------------------------------

TRUE_THETA3 <- c(delta = 0.5, alpha1 = 2, alpha2 = 3, beta1 = 1, beta2 = 1.2)
N <- 800
set.seed(200)
Data3 <- rgamma_mixture(N, TRUE_THETA3[8], TRUE_THETA3[9], TRUE_THETA3[10], TRUE_THETA3[11], TRUE_THETA3[12])
INITIAL_THETA3 <- c(0.5, 1, 1, 0.5, 0.5)

cat("\n--- Cenário 3: Mistura sobreposta (N=800) ---\n")
ECM_Result3 <- ECM_Gamma(Data3, INITIAL_THETA3)
print(paste("Iterações para convergência:", ECM_Result3$iterations))
df_estimates3 <- data.frame(
  Parameter = c("delta", "alpha1", "alpha2", "beta1", "beta2"),
  True = TRUE_THETA3,
  Estimate = ECM_Result3$estimates
)
print(round(df_estimates3, 5))