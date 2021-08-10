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