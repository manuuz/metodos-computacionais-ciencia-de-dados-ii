ic_bootstrap_studentizado <- function(x, stat = mean, B = 1000, R = 200, alpha = 0.05) {
  n <- length(x)
  theta_chapeu <- stat(x)
  ep_chapeu <- sd(replicate(R, stat(sample(x, n, replace = TRUE))))
  
  t_boot <- numeric(B)
  for (b in 1:B) {
    x_estrela <- sample(x, n, replace = TRUE)
    theta_b <- stat(x_estrela)
    ep_b <- sd(replicate(R, stat(sample(x_estrela, n, replace = TRUE))))
    t_boot[b] <- (theta_b - theta_chapeu) / ep_b
  }
  
  t_quant <- quantile(t_boot, c(alpha / 2, 1 - alpha / 2))
  ic <- c(theta_chapeu - t_quant[2] * ep_chapeu,
          theta_chapeu - t_quant[1] * ep_chapeu)
  return(ic)
}

ic_bootstrap_bca <- function(x, stat = mean, B = 1000, alpha = 0.05) {
  n <- length(x)
  theta_chapeu <- stat(x)
  
  # bootstrap
  thetas <- replicate(B, stat(sample(x, n, replace = TRUE)))
  
  # jackknife
  jack <- sapply(1:n, function(i) stat(x[-i]))
  theta_barra <- mean(jack)
  
  # parâmetros BCa
  z0 <- qnorm(mean(thetas < theta_chapeu))
  a <- sum((theta_barra - jack)^3) / (6 * (sum((theta_barra - jack)^2))^(3/2))
  
  z_alfa <- qnorm(c(alpha / 2, 1 - alpha / 2))
  alfa <- pnorm(z0 + (z0 + z_alfa) / (1 - a * (z0 + z_alfa)))
  
  ic <- quantile(thetas, probs = alfa)
  
  return(ic)
}

library(boot)
set.seed(123)
x <- rnorm(30, mean = 5, sd = 2)

# necessário para "stud"

metricas_boot_student <- function(data, idx) {
  d <- data[idx]
  m <- mean(d)
  se <- sd(d) / sqrt(length(d))
  return(c(mean = m, se = se))
}

# retorna média e erro padrão
boot_obj <- boot(x, statistic = metricas_boot_student, R = 1000, sim = "ordinary", stype = "i")


ic_student <- ic_bootstrap_studentizado(x, stat = mean, B = 500, R = 100)
ic_bca <- ic_bootstrap_bca(x, stat = mean, B = 1000)

# funções do pacote
boot_student <- boot.ci(boot_obj, type = "stud")
boot_bca <- boot.ci(boot_obj, type = "bca")

resultado <- data.frame(
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

print(resultado)