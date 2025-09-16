library(pracma)
n <- 10000

# exp{e^x} dx

f <- function(x){
  exp(exp(x))
}

x <- runif(n)
aprox <- mean(f(x))
print(aprox)

print(integrate(f = function(x){exp(exp(x))}, lower=0, upper=1))

# exp{x+x^2}

f <- function(x){
  exp(x-x^2)
}

u <- runif(n, -2, 2)
aprox <- mean(f(u)) * (2 - (-2))
print(aprox)

print(integrate(f = function(x){exp(x-x^2)}, lower=-2, upper=2))

# x(1+x^2)^-2

g <- function(u){
  (u*(1-u))/(2*u^2 - 2*u + 1)^2
}

u <- runif(n)
aprox <- mean(g(u))
print(aprox)

print(integrate(f = function(x){x*(1+x^2)^(-2)}, lower=0, upper=Inf))

# exp{-(x+y)}

g <- function(u, v) {
  return(u^v * (-log(u)))
}

u <- runif(n)
v <- runif(n)

aprox <- mean(g(u, v))
print(aprox)

f <- function(x, y) {
  exp(-(x + y))
}

print(integral2(f, 0, 100, 0, function(x){x}))
