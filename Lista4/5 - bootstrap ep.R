M <- 1000
B <- 2000
tamanhos <- c(30, 100, 500) 

r <- 4 
p <- 0.75 

var_x <- r * (1 - p) / (p^2) 

# calcula o erro padrão verdadeiro de X barra
calcula_ep <- function(n) {
  sqrt(var_x / n)
}

ep <- list()

simulacao <- function(tamanho_amostral, M, B, r, p) {
  
  ep_np <- numeric(M) # não paramétrico
  ep_p <- numeric(M) # paramétrico
  
  ep_verdadeiro <- calcula_ep(tamanho_amostral)
  
  cat(sprintf("Rodando para n = %d (SE Real: %.4f)...\n", tamanho_amostral, ep_verdadeiro))
  
  for (i in 1:M) {
    # geração da amostra real
    amostra <- rnbinom(tamanho_amostral, size = r, prob = p)
    
    # não paramétrico
    media_np <- replicate(B, {
      mean(sample(amostra, tamanho_amostral, replace = TRUE))
    })
    ep_np[i] <- sd(media_np)
    
    # paramétrico (mal especificado)
    lambda_chapeu <- mean(amostra)
    
    media_p <- replicate(B, {
      mean(rpois(tamanho_amostral, lambda = lambda_chapeu))
    })
    ep_p[i] <- sd(media_p)
  }
  
  return(data.frame(
    n = tamanho_amostral,
    EP_Verdadeiro = ep_verdadeiro,
    EP_NP = ep_np,
    EP_P = ep_p
  ))
}

set.seed(123456789)
for (n in tamanhos) {
  ep[[as.character(n)]] <- simulacao(n, M, B, r, p)
}

dados_plot <- data.frame()
for (n_str in names(ep)) {
  res <- ep[[n_str]]
  n <- res$n[1]  # apenas o valor de n (é constante dentro de cada data.frame)
  ep_real <- res$EP_Verdadeiro[1]
  
  df_np <- data.frame(
    n = factor(n),
    Metodo = "Não-Paramétrico",
    EP_Estimado = res$EP_NP,
    EP_Verdadeiro = ep_real
  )
  
  df_p <- data.frame(
    n = factor(n),
    Metodo = "Paramétrico (ajuste ruim)",
    EP_Estimado = res$EP_P,
    EP_Verdadeiro = ep_real
  )
  
  dados_plot <- rbind(dados_plot, df_np, df_p)
}

ggplot(dados_plot, aes(x = EP_Estimado, fill = Metodo)) +
  geom_density(alpha = 0.6) +
  facet_wrap(~ n, scales = "free_x", labeller = label_both) +
  geom_vline(aes(xintercept = EP_Verdadeiro),
             linetype = "dashed", color = "black") +
  scale_fill_manual(values = c(
    "Não-Paramétrico" = "blue",
    "Paramétrico (ajuste ruim)" = "red"
  )) +
  labs(
    title = "Comparação dos Estimadores de Erro Padrão por Bootstrap",
    subtitle = "Linha tracejada: erro padrão real",
    x = "EP Estimado",
    y = "Densidade"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")