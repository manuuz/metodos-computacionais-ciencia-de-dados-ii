# parâmetros
set.seed(123456789)
k <- 100000
g <- function(x) exp(x)

# Monte Carlo simples
U_mc <- runif(k)
X_mc <- g(U_mc)
tc_mc <- mean(X_mc) # tc = theta chapéu
var_tc_mc <- var(X_mc) / k
ep_tc_mc <- sqrt(var_tc_mc)

# variável antitética
U_ant <- 1 - U_mc
Y_ant <- (g(U_mc) + g(U_ant)) / 2
tc_ant <- mean(Y_ant)
var_tc_ant <- var(Y_ant) / k
ep_tc_ant <- sqrt(var_tc_ant)

reducao <- ((var_tc_mc - var_tc_ant) / var_tc_mc) * 100 # em percentual
valor_exato <- exp(1) - 1

resultado <- data.frame(
  Método = c("Monte Carlo", "Antitética"),
  Estimativa = c(tc_mc, tc_ant),
  Variância = c(var_tc_mc, var_tc_ant),
  ErroPadrão = c(ep_tc_mc, ep_tc_ant)
)

kable(resultado)
cat("\nValor exato:", round(valor_exato, 7),
    "\nRedução percentual na variância:", round(reducao, 4), "%\n")