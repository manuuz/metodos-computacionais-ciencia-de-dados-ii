# solução numérica para alfa
resolver_alfa <- function(C_k, alfa_inicial = 1) {
  fun <- function(a) digamma(a) - C_k
  resultado <- tryCatch({
    uniroot(fun, interval = c(1e-8, 1e6), tol = 1e-8)$root
  }, error = function(e) {
    # newton-raphson
    a <- max(alfa_inicial, 1e-6)
    for (i in 1:100) {
      g <- digamma(a) - C_k
      tg <- trigamma(a)
      if (is.na(g) || is.na(tg) || tg <= 0) break
      a_novo <- a - g / tg
      if (a_novo <= 0) a_novo <- a / 2
      if (abs(a_novo - a) < 1e-8) { a <- a_novo; break }
      a <- a_novo
    }
    a
  })
  as.numeric(resultado)
}


ecm_gama <- function(X, parametros_iniciais, max_iter = 500, tol = 1e-6) {
  
  # parâmetros iniciais
  theta <- as.numeric(parametros_iniciais[1:5])
  names(theta) <- c("delta", "alfa1", "alfa2", "beta1", "beta2")
  n <- length(X)
  historico <- matrix(NA, nrow = max_iter, ncol = 5)
  colnames(historico) <- names(theta)
  
  for (t in 1:max_iter) {
    delta_t  <- theta["delta"]
    alfa1_t <- theta["alfa1"]
    alfa2_t <- theta["alfa2"]
    beta1_t  <- theta["beta1"]
    beta2_t  <- theta["beta2"]
    
    # passo E
    f1 <- dgamma(X, shape = alfa1_t, rate = beta1_t)
    f2 <- dgamma(X, shape = alfa2_t, rate = beta2_t)
    # p_i = P(Z_i = 1 | x_i, theta^(t)), onde Z_i é a variável latente com distribuição Bernoulli
    p_i <- (delta_t * f1) / (delta_t * f1 + (1 - delta_t) * f2)
    # para garantir estabilidade numérica
    p_i <- pmin(pmax(p_i, 1e-10), 1 - 1e-10)
    
    # passo CM1: maximização condicional em relação a delta
    delta_novo <- mean(p_i)
    
    # passo CM2: maximização condicional em relação a alfa1 e beta1
    somatorio_p <- sum(p_i)
    somatorio_x_p <- sum(p_i * X)
    if (somatorio_p > 0) {
      beta1_novo <- alfa1_t * somatorio_p / somatorio_x_p
      C1 <- (sum(p_i * log(X)) / somatorio_p) + log(beta1_novo)
      alfa1_novo <- resolver_alfa(C1, alfa1_t)
    } else {
      beta1_novo <- beta1_t
      alfa1_novo <- alfa1_t
    }
    
    # passo CM3: maximização condicional em relação a alfa2 e beta2
    somatorio_1mp <- n - somatorio_p
    somatorio_x_1mp <- sum(X) - somatorio_x_p
    if (somatorio_1mp > 0) {
      beta2_novo <- alfa2_t * somatorio_1mp / somatorio_x_1mp
      C2 <- (sum((1 - p_i) * log(X)) / somatorio_1mp) + log(beta2_novo)
      alfa2_novo <- resolver_alfa(C2, alfa2_t)
    } else {
      beta2_novo <- beta2_t
      alfa2_novo <- alfa2_t
    }
    
    # atualização dos parâmetros
    theta_novo <- c(delta_novo, alfa1_novo, alfa2_novo, beta1_novo, beta2_novo)
    names(theta_novo) <- names(theta)
    
    # critério de convergência
    crit <- max(abs(theta_novo - theta) / pmax(abs(theta), 1e-8))
    theta <- theta_novo
    historico[t, ] <- theta
    
    if (crit < tol) break
    
  }
  
  list(estimativas = theta, iteracoes = t, historico = historico[1:t, ])
}

# gerando os dados
mistura_gama <- function(n, delta, alfa1, alfa2, beta1, beta2) {
  Z <- rbinom(n, 1, delta)
  rgamma(n, shape = ifelse(Z == 1, alfa1, alfa2),
         rate = ifelse(Z == 1, beta1, beta2))
}

set.seed(123456789)

# cenário 1: mistura bem separada
parametros1 <- c(delta = 0.6, alfa1 = 5, alfa2 = 2, beta1 = 1, beta2 = 0.2)
dados1 <- mistura_gama(500, parametros1["delta"], parametros1["alfa1"],
                       parametros1["alfa2"], parametros1["beta1"], parametros1["beta2"])
resultado1 <- ecm_gama(dados1, c(0.5, 1, 1, 0.5, 0.5))

# cenário 2: mistura com maior influência da segunda distribuição
parametros2 <- c(delta = 0.2, alfa1 = 10, alfa2 = 5, beta1 = 2, beta2 = 0.5)
dados2 <- mistura_gama(1000, parametros2["delta"], parametros2["alfa1"],
                       parametros2["alfa2"], parametros2["beta1"], parametros2["beta2"])
resultado2 <- ecm_gama(dados2, c(0.5, 1, 1, 0.5, 0.5))

# cenário 3: mistura com densidades sobrepostas

parametros3 <- c(delta = 0.5, alfa1 = 2, alfa2 = 3, beta1 = 1, beta2 = 1.2)
dados3 <- mistura_gama(800, parametros3["delta"], parametros3["alfa1"],
                       parametros3["alfa2"], parametros3["beta1"], parametros3["beta2"])
resultado3 <- ecm_gama(dados3, c(0.5, 1, 1, 0.5, 0.5))

resultados <- data.frame(
  Cenário = c("Cenário 1: Separada", "Cenário 2: Desigual", "Cenário 3: Sobreposta"),
  Iterações = c(resultado1$iteracoes, resultado2$iteracoes, resultado3$iteracoes),
  round(rbind(resultado1$estimativas, resultado2$estimativas, resultado3$estimativas), 4)
)
kable(resultados)

plot_componentes <- function(dados, parametros_verdadeiros, parametros_estimados, titulo) {
  x <- seq(0, max(dados), length.out = 500)
  
  # densidades teóricas e estimadas
  df <- data.frame(
    x = x,
    f1_verdadeira = parametros_verdadeiros["delta"] * dgamma(x, parametros_verdadeiros["alfa1"], parametros_verdadeiros["beta1"]),
    f2_verdadeira = (1 - parametros_verdadeiros["delta"]) * dgamma(x, parametros_verdadeiros["alfa2"], parametros_verdadeiros["beta2"]),
    f1_estimada = parametros_estimados["delta"] * dgamma(x, parametros_estimados["alfa1"], parametros_estimados["beta1"]),
    f2_estimada = (1 - parametros_estimados["delta"]) * dgamma(x, parametros_estimados["alfa2"], parametros_estimados["beta2"])
  )
  
  ggplot() +
    geom_histogram(aes(x = dados, y = after_stat(density)),
                   bins = 40, fill = "grey85", color = "white", alpha = 0.6) +
    
    geom_line(data = df, aes(x = x, y = f1_verdadeira, color = "C1 Verdadeira"),
              linetype = "dashed", linewidth = 1) +
    geom_line(data = df, aes(x = x, y = f2_verdadeira, color = "C2 Verdadeira"),
              linetype = "dashed", linewidth = 1) +
    geom_line(data = df, aes(x = x, y = f1_estimada, color = "C1 Estimada"),
              linewidth = 1.1) +
    geom_line(data = df, aes(x = x, y = f2_estimada, color = "C2 Estimada"),
              linewidth = 1.1) +
    
    scale_color_manual(
      name = "Curvas de densidade",
      values = c("C1 Verdadeira" = "blue",
                 "C2 Verdadeira" = "red",
                 "C1 Estimada" = "darkblue",
                 "C2 Estimada" = "darkred")
    ) +
    
    theme_minimal() +
    labs(title = titulo, x = "x", y = "Densidade") +
    theme(plot.title = element_text(size = 12, face = "bold"),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 9))
}


pairs = c(mfrow=c(1,3))
plot_componentes(dados1, parametros1, resultado1$estimativas, "Cenário 1")
plot_componentes(dados2, parametros2, resultado2$estimativas, "Cenário 2")
plot_componentes(dados3, parametros3, resultado3$estimativas, "Cenário 3")