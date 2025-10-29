# densidade gama (parametrização forma/taxa)
dgamma_mix <- function(x, alpha, beta) dgamma(x, shape = alpha, rate = beta)

# solução numérica para alfa via densidade da gama inversa
solve_alpha <- function(C_k, alpha_init = 1) {
  fun <- function(a) digamma(a) - C_k
  lower <- 1e-8; upper <- 1e6
  res <- tryCatch({
    uniroot(fun, interval = c(lower, upper), tol = 1e-8)$root
  }, error = function(e) {
    a <- max(alpha_init, 1e-6)
    for (i in 1:100) {
      g <- digamma(a) - C_k
      tg <- trigamma(a)
      if (is.na(g) || is.na(tg) || tg <= 0) break
      a_new <- a - g / tg
      if (a_new <= 0) a_new <- a / 2
      if (abs(a_new - a) < 1e-8) { a <- a_new; break }
      a <- a_new
    }
    a
  })
  as.numeric(res)
}

ECM_Gamma <- function(X, initial_params, max_iter = 500, tol = 1e-6) {
  theta <- as.numeric(initial_params[1:5])
  names(theta) <- c("delta", "alpha1", "alpha2", "beta1", "beta2")
  n <- length(X)
  history <- matrix(NA, nrow = max_iter, ncol = 5)
  colnames(history) <- names(theta)
  crit <- Inf
  
  for (t in 1:max_iter) {
    delta_t  <- theta["delta"]
    alpha1_t <- theta["alpha1"]
    alpha2_t <- theta["alpha2"]
    beta1_t  <- theta["beta1"]
    beta2_t  <- theta["beta2"]
    
    f1 <- dgamma_mix(X, alpha1_t, beta1_t)
    f2 <- dgamma_mix(X, alpha2_t, beta2_t)
    numerator <- delta_t * f1
    denominator <- numerator + (1 - delta_t) * f2
    p_i <- ifelse(denominator <= 0, 0, numerator / denominator)
    p_i <- pmin(pmax(p_i, 0), 1)
    
    sum_p <- sum(p_i); sum_1_minus_p <- n - sum_p
    sum_x_p <- sum(X * p_i); sum_x_1_minus_p <- sum(X) - sum_x_p
    
    delta_new <- sum_p / n
    beta1_new <- if (sum_p > 0) alpha1_t * sum_p / sum_x_p else beta1_t
    beta2_new <- if (sum_1_minus_p > 0) alpha2_t * sum_1_minus_p / sum_x_1_minus_p else beta2_t
    
    alpha1_new <- if (sum_p > 0) solve_alpha((sum(p_i * log(X)) / sum_p) + log(beta1_new), alpha1_t) else alpha1_t
    alpha2_new <- if (sum_1_minus_p > 0) solve_alpha((sum((1 - p_i) * log(X)) / sum_1_minus_p) + log(beta2_new), alpha2_t) else alpha2_t
    
    theta_new <- c(delta_new, alpha1_new, alpha2_new, beta1_new, beta2_new)
    names(theta_new) <- names(theta)
    
    crit <- max(abs(theta_new - theta) / pmax(abs(theta), 1e-8))
    theta <- theta_new
    history[t, ] <- theta
    if (crit < tol) break
  }
  list(estimates = theta, iterations = t)
}

# gerando dados da mistura
rgamma_mixture <- function(n, delta, alpha1, alpha2, beta1, beta2) {
  Z <- rbinom(n, 1, delta)
  rgamma(n, shape = ifelse(Z == 1, alpha1, alpha2), rate = ifelse(Z == 1, beta1, beta2))
}

set.seed(123456789)

# cenário 1: mistura bem separada (parâmetros de escala e forma diferentes)
TRUE_THETA1 <- c(delta = 0.6, alpha1 = 5, alpha2 = 2, beta1 = 1, beta2 = 0.2)
Data1 <- rgamma_mixture(500, TRUE_THETA1["delta"], TRUE_THETA1["alpha1"], TRUE_THETA1["alpha2"], TRUE_THETA1["beta1"], TRUE_THETA1["beta2"])
res1 <- ECM_Gamma(Data1, c(0.5, 1, 1, 0.5, 0.5))

# cenário 2: mistura com tamanhos desiguais (maior peso no grupo 2)
TRUE_THETA2 <- c(delta = 0.2, alpha1 = 10, alpha2 = 5, beta1 = 2, beta2 = 0.5)
Data2 <- rgamma_mixture(1000, TRUE_THETA2["delta"], TRUE_THETA2["alpha1"], TRUE_THETA2["alpha2"], TRUE_THETA2["beta1"], TRUE_THETA2["beta2"])
res2 <- ECM_Gamma(Data2, c(0.5, 1, 1, 0.5, 0.5))

# cenário 3: mistura com densidades mais sobrepostas
TRUE_THETA3 <- c(delta = 0.5, alpha1 = 2, alpha2 = 3, beta1 = 1, beta2 = 1.2)
Data3 <- rgamma_mixture(800, TRUE_THETA3["delta"], TRUE_THETA3["alpha1"], TRUE_THETA3["alpha2"], TRUE_THETA3["beta1"], TRUE_THETA3["beta2"])
res3 <- ECM_Gamma(Data3, c(0.5, 1, 1, 0.5, 0.5))

resultados <- data.frame(
  Cenário = c("Cenário 1: Separada", "Cenário 2: Desigual", "Cenário 3: Sobreposta"),
  Iterações = c(res1$iterations, res2$iterations, res3$iterations),
  round(rbind(res1$estimates, res2$estimates, res3$estimates), 4)
)

kable(resultados)