source("simulation functions.R")

dat <- read.csv("simulated.csv")
dat <- dat[!is.na(dat$N_est), c("N_est", "N_pres", "N_det")]

# Simulate the study site for three years, with N expected individuals

transect_list <- read.traps(file = "transects.txt", detector = "transect")
x <- transect_list$x
y <- transect_list$y

x0 <- x[seq(1, 22, 2)]
y0 <- y[seq(1, 22, 2)]
x1 <- x[seq(2, 22, 2)]
y1 <- y[seq(2, 22, 2)]

pars <- list(lambda0 = 1, sigma = 30)
mask <- make.mask(traps = transect_list, buffer = 0, spacing = 10)
times <- c(36500, 36530, 36560)

inits <- function(){
  list(beta0 = rnorm(1),
       beta1 = rnorm(1),
       sig2 = runif(1))
}

monitor <- c("beta0", "beta1", "sig2", "predict", "corrected")

nc <- 3
nb <- 5000
ni <- 10000 + nb
nt <- 1

B <- 5
expected <- numeric(B)
present <- numeric(B)
detectable <- numeric(B)
uncorrected <- numeric(B)
corrected <- numeric(B)
median <- numeric(B)
lcl <- numeric(B)
ucl <- numeric(B)

for(i in 1:B){
  N_exp <- 30
  pop <- population(30, 365 * 101, MEAN_LIFE)
  capthist <- survey_sim(x0, y0, x1, y1, pop, times = times, HHN, pars, t0=36500-30)
  # mod <- secr.fit(capthist = capthist$secr_capthist, mask = mask)
  # 
  # jagsdat <- list(N = nrow(dat),
  #                 present = dat$N_pres,
  #                 estimated = dat$N_est,
  #                 point_est = region.N(mod)[2, 1],
  #                 std_error = region.N(mod)[2, 2])
  # 
  # out <- jags(data = jagsdat,
  #             inits = inits,
  #             parameters.to.save = monitor,
  #             model.file = "model.txt",
  #             n.chains = nc,
  #             n.iter = ni,
  #             n.burnin = nb,
  #             n.thin = nt)
  # 
  expected[i] <- N_exp
  present[i] <- pres_det(pop, 36500 - 30)$pres
  detectable[i] <- pres_det(pop, 36500 - 30)$det
  # uncorrected[i] <- jagsdat$point_est
  # corrected[i] <- out$mean$corrected
  # median[i] <- out$q50$corrected
  # lcl[i] <- out$q2.5$corrected
  # ucl[i] <- out$q97.5$corrected
  uncorrected[i] <- survreg(make_surv_object(capthist$survhist, 36500-30, times) ~ 1) %>% predict()[1]
}

plot(present, corrected, ylim = range(lcl, ucl))
segments(present, lcl, present, ucl)
abline(0, 1)

beep()
