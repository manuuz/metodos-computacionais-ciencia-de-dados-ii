
require(microbenchmark)
require(ggplot2)

rpoisson_aux_livro <- function(lambda){
  u <- runif(1, 0, 1)
  i <- 0; pr <- exp(-lambda); Fx = pr
  while(u >= Fx){
    pr <- pr*lambda/(i+1)
    Fx <- Fx + pr
    i <- i+1
  }
  return(i)
}

rpoisson_livro <- function(n, lambda){
  replicate(n, expr = rpoisson_aux_livro(lambda), simplify = TRUE)
}


rpoisson_aux_knuth <- function(lambda){
  L <- exp(-lambda); k <- 0; p <- 1
  
  while(p > L){
    k <- k+1
    u <- runif(1, 0, 1)
    p <- p*u
  }
  
  return(k-1)
}


rpoisson_knuth <- function(n, lambda){
  replicate(n, expr = rpoisson_aux_knuth(lambda), simplify = TRUE)
}

n <- 1000 
lambda <- 5 

mb_results <- microbenchmark(
  Livro_Referência = rpoisson_livro(n, lambda),
  Otimizado_Knuth = rpoisson_knuth(n, lambda),
  Nativo_R = rpois(n, lambda),
  times = 100L 
)

print(mb_results)

autoplot(mb_results) +
  labs(title = paste("Desempenho para lambda =", lambda)) +
  theme_minimal()

