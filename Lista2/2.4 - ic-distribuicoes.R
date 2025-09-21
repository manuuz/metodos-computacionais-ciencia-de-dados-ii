set.seed(2023032088)
m <- 20000
alfa <- 0.05
confianca_nominal <- 1 - alfa
n <- c(20, 30, 100)
resultados <- list()
i <- 1

estimar <- function(coberturas, m, dist, n) {
  y_chapeu <- mean(coberturas)
  ep_estimado <- sqrt(y_chapeu * (1 - y_chapeu) / m)
  
  return(data.frame(
    Distribuicao = dist,
    n = n,
    Cobertura_Empirica = y_chapeu,
    Erro_Padrao = ep_estimado)
    )
}

# a
mu_chisq <- 2

for (n in n) {
  coberturas_chisq <- replicate(m, expr = {
    x <- rchisq(n, df = 2)
    teste_t <- t.test(x, conf.level = confianca_nominal)
    confianca <- teste_t$conf.int
    confianca[1] <= mu_chisq && mu_chisq <= confianca[2]
    }
  )

  resultados[[i]] <- estimar(coberturas_chisq, m, "Chi^2_2 (mu=2)", n)
  i <- i + 1
}

# b
mu_gamma <- 2 

for (n in n) {
  coberturas_gamma <- replicate(m, expr = {
    x <- rgamma(n, shape = 2, rate = 1) # Gama(2,1)
    teste_t <- t.test(x, conf.level = confianca_nominal)
    confianca <- teste_t$conf.int
    confianca[1] <= mu_gamma && mu_gamma <= confianca[2]
    }
  )
  resultados[[i]] <- estimar(coberturas_gamma, m, "Gamma(2, 1) (mu=2)", n)
  i <- i + 1
}

tabela_final <- do.call(rbind, resultados)
tabela_final$Cobertura_Empirica <- round(tabela_final$Cobertura_Empirica, 4)
tabela_final$Erro_Padrao <- round(tabela_final$Erro_Padrao, 5)
kable(tabela_final, caption = "IC de 95% para distribuições assimétricas")