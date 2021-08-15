### Test Data ###
library(survival)
library(npsurv)
data(leukemia)

# Create a data frame with data for testing
t <- leukemia[, 1]
d <- as.integer(leukemia[, 2] < Inf)
tmt <- c(rep(1, 21), rep(0, 21))

dat <- data.frame(t=t, d=d, tmt=tmt)

### Exponential ###
# Negative log-likelihood function for the exponential response
exp.nll <- function(theta, data){
  # Select which observations were observed
  obs <- which(data$d == 1)
  
  # Calculate parameter
  lambda <- 1 / exp(theta[1] + theta[2] * data$tmt)
  
  # Find likelihood of observed times
  ll <- sum(dexp(data$t[obs], lambda[obs], log=TRUE))
  # Add likelihood of censored times
  ll <- ll + sum(pexp(data$t[-obs], lambda[-obs], lower.tail = FALSE, log.p = TRUE))
  
  return(-ll)
}

# Fit model with survreg and compare to manual optimisation
exp.fit <- survreg(Surv(t, d) ~ tmt, dist = "exponential")
summary(exp.fit)

exp.mle <- optim(rnorm(2), exp.nll, data = dat)
exp.mle$par

### Weibull ###
# Weibull density and cumulative density functions using alternate parameterisation
dweib <- function(t, shape, scale, log = FALSE){
  rescale <- (1 / scale) ^ (1 / shape)
  dweibull(t, shape, rescale, log)
}

pweib <- function(t, shape, scale, lower.tail = TRUE, log.p = FALSE){
  rescale <- (1 / scale) ^ (1 / shape)
  pweibull(t, shape, rescale, lower.tail, log.p)
}

# Negative log-likelihood functions for the Weibull response
weib.nll <- function(theta, data){
  # Select which observations were observed
  obs <- which(d == 1)
  
  # Shape parameter
  p <- exp(theta[1])
  # Scale parameter
  lambda <- exp(theta[2] + theta[3] * data$tmt) ^ (- p)
  
  # Find likelihood of observed times
  ll <- sum(dweib(data$t[obs], p, lambda[obs], log = TRUE))
  # Add likelihood of censored times
  ll <- ll + sum(pweib(data$t[-obs], p, lambda[-obs], lower.tail = FALSE, log.p = TRUE))
  return(-ll)
}

# Fit model with survreg and compare to manual optimisation
weib.fit <- survreg(Surv(t, d) ~ tmt, dist = "weibull")
summary(weib.fit)

weib.mle <- optim(rnorm(3), weib.nll, data = dat)
weib.mle$par

### Log-logistic ###
# Density and cumulative density functions for the log-logistic distribution
dloglog <- function(x, shape, scale, log = FALSE){
  density <- (scale * shape * (scale * x) ** (shape - 1)) / ((1 + (scale * x) ** shape) ** 2)
  
  if(log){
    log(density)
  }
  else{
    density
  }
}

ploglog <- function(q, shape, scale, lower.tail = TRUE, log.p = FALSE){
  prob <- ((scale * q) ** shape) / (1 + ((scale * q) ** shape))
  
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

# Negative log-likelihood function for the log logistic response
loglog.nll <- function(theta, data){
  # Select which observations were observed
  obs <- which(data$d == 1)
  
  # Calculate parameters
  shape <- exp(theta[1])
  scale <- exp(theta[2] + theta[3] * data$tmt) ^ -1
  
  # Find likelihood of observed times
  ll <- sum(dloglog(data$t[obs], shape, scale[obs], log = TRUE))
  # Add likelihood of censored times
  ll <- ll + sum(ploglog(data$t[-obs], shape, scale[-obs], lower.tail = FALSE, log.p = TRUE))
  
  return(-ll)
}

# Fit model with survreg and compare to manual optimisation
loglog.fit <- survreg(Surv(t, d) ~ tmt, dist = "loglogistic")
summary(loglog.fit)

loglog.mle <- optim(rnorm(3), loglog.nll, data = dat)
loglog.mle$par
