ic_bootstrap_studentizado <- function(x, stat = mean, B = 1000, R = 200, alpha = 0.05) {
  n <- length(x)
  theta_hat <- stat(x)
  se_hat <- sd(replicate(R, stat(sample(x, n, replace = TRUE))))
  
  t_boot <- numeric(B)
  for (b in 1:B) {
    x_star <- sample(x, n, replace = TRUE)
    theta_b <- stat(x_star)
    se_b <- sd(replicate(R, stat(sample(x_star, n, replace = TRUE))))
    t_boot[b] <- (theta_b - theta_hat) / se_b
  }
  
  t_quant <- quantile(t_boot, c(alpha / 2, 1 - alpha / 2))
  ci <- c(theta_hat - t_quant[2] * se_hat,
          theta_hat - t_quant[1] * se_hat)
  return(ci)
}

ic_bootstrap_bca <- function(x, stat = mean, B = 1000, alpha = 0.05) {
  n <- length(x)
  theta_hat <- stat(x)
  
  # Bootstrap
  
  thetas <- replicate(B, stat(sample(x, n, replace = TRUE)))
  
  # Jackknife
  
  jack_vals <- sapply(1:n, function(i) stat(x[-i]))
  theta_bar <- mean(jack_vals)
  
  # Parâmetros BCa
  
  z0 <- qnorm(mean(thetas < theta_hat))
  a <- sum((theta_bar - jack_vals)^3) / (6 * (sum((theta_bar - jack_vals)^2))^(3/2))
  
  z_alpha <- qnorm(c(alpha / 2, 1 - alpha / 2))
  alpha_adj <- pnorm(z0 + (z0 + z_alpha) / (1 - a * (z0 + z_alpha)))
  
  ci <- quantile(thetas, probs = alpha_adj)
  return(ci)
}

library(boot)

#Dados de exemplo

set.seed(123)
x <- rnorm(30, mean = 5, sd = 2)

#Função de estatística que retorna valor e erro padrão (necessário para "stud")

mean_boot_student <- function(data, idx) {
  d <- data[idx]
  m <- mean(d)
  se <- sd(d) / sqrt(length(d))
  return(c(mean = m, se = se))
}

#Bootstrap com a função que retorna mean + se

boot_obj <- boot(x, statistic = mean_boot_student, R = 1000, sim = "ordinary", stype = "i")

#Intervalos via funções próprias

ic_student <- ic_bootstrap_studentizado(x, stat = mean, B = 500, R = 100)
ic_bca <- ic_bootstrap_bca(x, stat = mean, B = 1000)

#Intervalos com boot.ci()

boot_student <- boot.ci(boot_obj, type = "stud")
boot_bca <- boot.ci(boot_obj, type = "bca")

df_ic <- data.frame(
  Metodo = c("Studentizado (manual)",
             "BCa (manual)",
             "Studentizado (pacote)",
             "BCa (pacote)"),
  IC_Inferior = c(
    round(ic_student[1], 4),
    round(ic_bca[1], 4),
    round(as.numeric(boot_student$student[,4]), 4),
    round(boot_bca$bca[,4], 4)
  ),
  IC_Superior = c(
    round(ic_student[2], 4),
    round(ic_bca[2], 4),
    round(as.numeric(boot_student$student[,5]), 4),
    round(boot_bca$bca[,5], 4)
  )
)

print(df_ic)