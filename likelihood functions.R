### Exponential ###
nll <- function(theta, data){
  obs <- which(data$d == 1)
  
  lambda <- 1 / exp(theta[1] + theta[2] * data$tmt)
  ll <- sum(dexp(data$t[obs], lambda[obs], log=TRUE))
  ll <- ll + sum(pexp(data$t[-obs], lambda[-obs], lower.tail = FALSE, log.p = TRUE))
  
  return(-ll)
}

### Weibull ###
# Weibull density function using alternate parameterisation
dweib <- function(t, shape, scale, log = FALSE){
  rescale <- (1 / scale) ^ (1 / shape)
  dweibull(t, shape, rescale, log)
}

pweib <- function(t, shape, scale, lower.tail = TRUE, log.p = FALSE){
  rescale <- (1 / scale) ^ (1 / shape)
  pweibull(t, shape, rescale, lower.tail, log.p)
}

nll <- function(theta, data){
  obs <- which(d == 1)
  
  p <- exp(theta[1])
  lambda <- exp(theta[2] + theta[3] * data$tmt) ^ (- p)
  
  ll <- sum(dweib(data$t[obs], p, lambda[obs], log = TRUE))
  ll <- ll + sum(pweib(data$t[-obs], p, lambda[-obs], lower.tail = FALSE, log.p = TRUE))
  return(-ll)
}

### Log-logistic ###
dloglog <- function(x, shape, scale, log = FALSE){
  density <- (scale * shape * x ** (shape - 1)) / ((1 + scale * x ** shape) ** 2)
  
  if(log){
    log(density)
  }
  else{
    density
  }
}

ploglog <- function(q, shape, scale, lower.tail = TRUE, log.p = FALSE){
  prob <- 1 / (1 + (scale * q ** shape) ** -1)
  
  if(!lower.tail){
    prob <- 1 - prob
  }
  else{
  }
  
  if(log.p){
    log(prob)
  }
  else{
    prob
  }
}

nll <- function(theta, data){
  obs <- which(data$d == 1)
  
  shape <- exp(theta[1])
  scale <- exp(theta[2] + theta[3] * data$tmt) ^ (-shape)
  
  ll <- sum(dloglog(data$t[obs], shape, scale[obs], log = TRUE))
  ll <- ll + sum(ploglog(data$t[-obs], shape, scale[-obs], lower.tail = FALSE, log.p = TRUE))
  
  return(-ll)
}

library(survival)
library(npsurv)
data(leukemia)

t <- leukemia[, 1]
d <- as.integer(leukemia[, 2] < Inf)
tmt <- c(rep(1, 21), rep(0, 21))

dat <- data.frame(t=t, d=d, tmt=tmt)

mle <- optim(rnorm(3), nll, data = dat, hessian = TRUE)
mle$par

fit <- survreg(Surv(t, d) ~ tmt)
summary(fit)