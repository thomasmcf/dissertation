library(jagsUI)
library(MCMCvis)

dat <- read.csv("simulated.csv")
dat <- dat[!is.na(dat$N_est), c("N_est", "N_pres", "N_det")]

inits <- function(){
  list(beta0 = rnorm(1),
       beta1 = rnorm(1),
       sig2 = runif(1))
}

jagsdat <- list(N = nrow(dat),
                present = dat$N_pres,
                estimated = dat$N_est,
                point_est = region.N(mod)[2, 1],
                std_error = region.N(mod)[2, 2])

monitor <- c("beta0", "beta1", "sig2", "predict", "corrected")

nc <- 3
nb <- 5000
ni <- 10000 + nb
nt <- 1

out <- jags(data = jagsdat,
            inits = inits,
            parameters.to.save = monitor,
            model.file = "model.txt",
            n.chains = nc,
            n.iter = ni,
            n.burnin = nb,
            n.thin = nt)

with(jagsdat, plot(pre, det))
abline(out$q50$beta0, out$q50$beta1)

lcl <- apply(out$sims.list$predict, 2, quantile, probs = 0.025)
ucl <- apply(out$sims.list$predict, 2, quantile, probs = 0.975)

lines(1:200, lcl)
lines(1:200, ucl)
