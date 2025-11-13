teste_boot <- function(z, mu0, B = 1000) {
  n <- length(z)
  
  teste_t_obs <- t.test(z, mu = mu0)
  t_obs <- teste_t_obs$statistic
  
  # translação para ter média mu0.
  z_til <- z - mean(z) + mu0
  
  t_estrela <- numeric(B)
  
  for(b in 1:B){
    # amostra centrada em mu0
    z_estrela <- sample(z_til, size = n, replace = TRUE)
    t_estrela[b] <- t.test(z_estrela, mu = mu0)$statistic
  }
  
  p_valor <- ( sum(t_estrela >= abs(t_obs)) + sum(t_estrela <= -abs(t_obs)) ) / B
  
  return(p_valor)
}

poder_mc <- function(mu_vec, n, M, B, alfa) {
  
  num_mu <- length(mu_vec)
  poder_t <- numeric(num_mu)
  poder_boot <- numeric(num_mu)
  mu0 <- 2
  
  for (i in 1:num_mu) {
    mu_real <- mu_vec[i]
    # Em uma distribuicao Chi-quadrado, E[X] = nu. Entao, nu = mu_real.
    nu_real <- mu_real 
    
    rejeicoes_t <- 0
    rejeicoes_boot <- 0
    
    # Monte Carlo
    for (m in 1:M) {
      dados_mc <- rchisq(n, df = nu_real)
      
      p_valor_t <- t.test(dados_mc, mu = mu0, alternative = "two.sided")
      if (p_valor_t$p.value <= alfa) {
        rejeicoes_t <- rejeicoes_t + 1
      }
      
      p_valor_boot <- teste_boot(dados_mc, mu0 = mu0, B = B)
      if (p_valor_boot <= alfa) {
        rejeicoes_boot <- rejeicoes_boot + 1
      }
    }
    
    poder_t[i] <- rejeicoes_t / M
    poder_boot[i] <- rejeicoes_boot / M
  }
  
  return(data.frame(
    mu_real = mu_vec,
    Poder_T = poder_t,
    Poder_Bootstrap = poder_boot
  ))
}

mu_valores <- c(1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0)
M <- 500
B <- 1000
alfa <- 0.05

set.seed(123456789)

resultados1 <- poder_mc(
  mu_vec = mu_valores,
  n = 30,
  M = M,
  B = B,
  alfa = alfa
)

kable(resultados1)

df_plot1 <- data.frame(
  mu_real = resultados1$mu_real,
  Poder = c(resultados1$Poder_T, resultados1$Poder_Bootstrap),
  Metodo = factor(rep(c("Teste T (Violado)", "Bootstrap"), each = length(mu_valores)))
)

resultados2 <- poder_mc(
  mu_vec = mu_valores,
  n = 100,
  M = M,
  B = B,
  alfa = alfa
)

kable(resultados2)

par(mfrow=c(1,2))

df_plot2 <- data.frame(
  mu_real = resultados2$mu_real,
  Poder = c(resultados2$Poder_T, resultados2$Poder_Bootstrap),
  Metodo = factor(rep(c("Teste T (Violado)", "Bootstrap"), each = length(mu_valores)))
)

ggplot(df_plot1, aes(x = mu_real, y = Poder, color = Metodo)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_vline(xintercept = 2, linetype = "dashed", color = "black", alpha = 0.7) +
  geom_hline(yintercept = alfa, linetype = "dotted", color = "gray50") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Poder: Teste T vs. Bootstrap (n = 30)",
    subtitle = "População Qui-quadrado",
    x = expression(paste(mu)),
    y = "Poder"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggplot(df_plot2, aes(x = mu_real, y = Poder, color = Metodo)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_vline(xintercept = 2, linetype = "dashed", color = "black", alpha = 0.7) +
  geom_hline(yintercept = alfa, linetype = "dotted", color = "gray50") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Poder: Teste T vs. Bootstrap (n = 100)",
    subtitle = "População Qui-quadrado",
    x = expression(paste(mu)),
    y = "Poder"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")