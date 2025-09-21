# Definindo parâmetros globais
M <- 20000          # Número de réplicas Monte Carlo. (Grande para erro padrão baixo)
alfa <- 0.05       # Nível de significância
conf_nominal <- 1 - alfa # Nível de confiança nominal (0.95)
n_valores <- c(20, 30, 100) # Tamanhos amostrais a serem testados

# Vetor para armazenar os resultados
resultados <- list()
i <- 1

# Função auxiliar para calcular a cobertura e o erro padrão do MC
calcular_resultados <- function(coberturas, M, nome_dist, n) {
  Y_hat <- mean(coberturas)
  SE_hat <- sqrt(Y_hat * (1 - Y_hat) / M)
  
  return(data.frame(
    Distribuicao = nome_dist,
    n = n,
    Cobertura_Empirica = Y_hat,
    Erro_Padrao_MC = SE_hat
  ))
}

# ------------------------------------------------------------------------------
# a) Distribuição Qui-quadrado com 2 graus de liberdade (Chi^2_2)
# Média populacional verdadeira (mu) = df = 2.
# ------------------------------------------------------------------------------

mu_chisq <- 2

for (n in n_valores) {
  # Simulação Monte Carlo
  # Para cada réplica, geramos a amostra, calculamos o IC t-Student de 95%,
  # e verificamos se a verdadeira média (mu_chisq = 2) está dentro.
  
  coberturas_chisq <- replicate(M, expr = {
    # Geração de amostra da Chi^2_2 (função rchisq)
    x <- rchisq(n, df = 2)
    
    # Cálculo do IC t-Student (utilizando a função t.test do R)
    teste_t <- t.test(x, conf.level = conf_nominal)
    CI <- teste_t$conf.int
    
    # Verificação da cobertura (y^j = 1 se coberto, 0 caso contrário)
    CI[4] <= mu_chisq && mu_chisq <= CI[5]
  })
  
  # Armazenando resultados
  resultados[[i]] <- calcular_resultados(coberturas_chisq, M, 
                                         "Chi^2_2 (mu=2)", n)
  i <- i + 1
}

# ------------------------------------------------------------------------------
# b) Outra distribuição assimétrica (Distribuição Gama)
# Escolha: Gama(shape=2, rate=1). Média populacional verdadeira (mu) = 2/1 = 2.
# ------------------------------------------------------------------------------

mu_gamma <- 2 

for (n in n_valores) {
  # Simulação Monte Carlo
  
  coberturas_gamma <- replicate(M, expr = {
    # Geração de amostra da Gama(shape=2, rate=1) (função rgamma)
    x <- rgamma(n, shape = 2, rate = 1)
    
    # Cálculo do IC t-Student
    teste_t <- t.test(x, conf.level = conf_nominal)
    CI <- teste_t$conf.int
    
    # Verificação da cobertura
    CI[4] <= mu_gamma && mu_gamma <= CI[5]
  })
  
  # Armazenando resultados
  resultados[[i]] <- calcular_resultados(coberturas_gamma, M, 
                                         "Gamma(2, 1) (mu=2)", n)
  i <- i + 1
}

# Combinando e exibindo a tabela final
tabela_final <- do.call(rbind, resultados)

# Arredondando para melhor visualização
tabela_final$Cobertura_Empirica <- round(tabela_final$Cobertura_Empirica, 4)
tabela_final$Erro_Padrao_MC <- round(tabela_final$Erro_Padrao_MC, 5)

print("Resultados do Estudo Monte Carlo (M = 20000 réplicas)")
print(tabela_final)