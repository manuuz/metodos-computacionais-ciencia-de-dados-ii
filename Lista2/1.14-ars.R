install.packages("ars")
library(ars)

# Gamma(2, 1/10)
gamma <-function(x,shape,scale=1){
  (shape-1)*log(x)-x/scale
  }

gamma_dlog <- function(x,shape,scale=1) {
  (shape-1)/x-1/scale
}
n<- 10000
sample_g<-ars(n, gamma, gamma_dlog, x=4.5, m=1, lb=TRUE, xlb=0, shape=2, scale=0.1)

hist(sample_g,
     probability = TRUE,
     ylim = c(0, 4),
     main = "Distribuição Gamma",
     col = 'lightblue')

x_seq <- seq(0, max(sample_g), length.out = 10000)
lines(x_seq, dgamma(x_seq, shape=2, scale=0.1), col = "tomato2", lwd = 3)


# Beta(3.4, 1.1)
beta <- function(x,a,b){
  (a-1)*log(x)+(b-1)*log(1-x)
}

beta_dlog <-function(x,a,b){
  (a-1)/x-(b-1)/(1-x)
}

sample_b<-ars(n, beta, beta_dlog, x=c(0.3,0.6), m=2, lb=TRUE, xlb=0, ub=TRUE, xub=1, a=3.4, b=1.1)

hist(sample_b,
     probability = TRUE,
     ylim = c(0, 4),
     main = "Distribuição Beta",
     col = 'lightblue')

x_seq <- seq(0, max(sample_b), length.out = 10000)
lines(x_seq, dbeta(x_seq, 3.4, 1.1), col = "tomato2", lwd = 3)

