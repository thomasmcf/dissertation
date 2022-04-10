source("simulation functions.R")

transect_list <- read.traps(file = "transects.txt", detector = "transect")
x <- transect_list$x
y <- transect_list$y

x0 <- x[seq(1, 22, 2)]
y0 <- y[seq(1, 22, 2)]
x1 <- x[seq(2, 22, 2)]
y1 <- y[seq(2, 22, 2)]

pars <- list(lambda0 = 1, sigma = 1)
mask <- make.mask(traps = transect_list, buffer = 0, spacing = 10)
times <- c(36500, 36530, 36560)
t0 <- 36500 - 30

cl <- makeCluster(2)
registerDoParallel(cl)
getDoParWorkers()

# Simulate populations of different abundances - with surveys
set.seed(012345)
Sys.time()

N_expected <- 1:5
B <- length(N_expected)

output <- foreach(i = 1:B, .packages = c("tidyverse", "secr", "survival"), .combine = rbind) %dopar% {
  N_exp <- N_expected[i]
  pop <- prune(population(N_exp, 365 * 101), 36500 - 30)
  
  capthist <- survey_sim(x0, y0, x1, y1, pop, times = times, t0=36500-30)
  mod <- secr::secr.fit(capthist = capthist$capthist, mask = mask)
  
  data.frame(expected = N_exp,
             present = pres_det(pop, 36500)$pres,
             detectable = pres_det(pop, 36500)$det,
             estimated = region.N(mod)[2, 1]
  )
}

beep()
saveRDS(output, "exp_pres_det_est.rds")



# Simulate populations with different mean dropping life time
mean_lifes <- 7:372
B <- length(mean_lifes)

output <- foreach(i = 1:B, .packages = c("tidyverse", "secr", "survival"), .combine = rbind) %dopar% {
  mean_life <- mean_lifes[i]
  pop <- population(30, 365 * 101, mean_life = mean_life)
  
  detectable_times <- sapply(pop,
                             function(ind) max(ind$degrade_times) - min(ind$drop_times))
  
  mean_detectable_time <- mean(detectable_times)
  se <- sd(detectable_times) / sqrt(length(detectable_times))
  
  slope <- mean_detectable_time / (365 + 2)
  
  capthist <- survey_sim(x0, y0, x1, y1, pop, times = times, t0=36500-30)
  mod <- secr::secr.fit(capthist = capthist$capthist, mask = mask)
  
  data.frame(estimated = region.N(mod)[2, 1],
             present = pres_det(pop, 36500)$pres,
             detectable = pres_det(pop, 36500)$det,
             mean_life = mean_life,
             mean_detectable_time = mean_detectable_time,
             se = se,
             slope = slope
  )
}

saveRDS(output, "dropping_life.rds")



# Simulate populations with different dropping rates
dropping_rate <- seq(1/7, 2, length.out = 200)
B <- length(dropping_rate)

output <- foreach(i = 1:B, .packages = c("tidyverse", "secr", "survival"), .combine = rbind) %dopar% {
  drop_rate <- dropping_rate[i]
  pop <- population(30, 365 * 101, drop_rate = drop_rate)
  
  detectable_times <- sapply(pop,
                             function(ind) max(ind$degrade_times) - min(ind$drop_times))
  
  mean_detectable_time <- mean(detectable_times)
  se <- sd(detectable_times) / sqrt(length(detectable_times))
  
  slope <- mean_detectable_time / (365 + (1 / drop_rate))
  
  capthist <- survey_sim(x0, y0, x1, y1, pop, times = times, t0=36500-30)
  mod <- secr::secr.fit(capthist = capthist$capthist, mask = mask)
  
  data.frame(drop_rate = drop_rate,
             mean_detectable_time = mean_detectable_time,
             se = se,
             slope = slope,
             estimated = region.N(mod)[2, 1],
             present = pres_det(pop, 36500)$pres,
             detectable = pres_det(pop, 36500)$det
  )
}

saveRDS(output, "drop_rate.rds")



# Simulate populations with different turnover rates
mean_stays <- seq(365/2, 365 * 1.5, length.out = 365)
B <- length(mean_stays)

output <- foreach(i = 1:B, .packages = c("tidyverse", "secr", "survival"), .combine = rbind) %dopar% {
  mean_stay <- mean_stays[i]
  pop <- population(30, 365 * 101, mean_stay = mean_stay)
  
  detectable_times <- sapply(pop,
                             function(ind) max(ind$degrade_times) - min(ind$drop_times))
  
  mean_detectable_time <- mean(detectable_times)
  se <- sd(detectable_times) / sqrt(length(detectable_times))
  
  slope <- mean_detectable_time / (mean_stay + 2)
  
  capthist <- survey_sim(x0, y0, x1, y1, pop, times = times, t0=36500-30)
  mod <- secr::secr.fit(capthist = capthist$capthist, mask = mask)
  
  data.frame(mean_stay = mean_stay,
             mean_detectable_time = mean_detectable_time,
             se = se,
             slope = slope,
             estimated = region.N(mod)[2, 1],
             present = pres_det(pop, 36500)$pres,
             detectable = pres_det(pop, 36500)$det
  )
}

saveRDS(output, "mean_stay.rds")

B <- 200
point_ests_3 <- foreach(i = 1:B, .packages = c("tidyverse", "secr", "survival"), .combine = c) %dopar% {
  pop <- prune(population(30, 36500), t0)
  
  survhist <- survey_sim(x0, y0, x1, y1, pop, times, t0)$survhist
  coef(survreg(survhist ~ 1, dist = "exponential"))
}

times <- 36500 + 0:11 * 30
point_ests_12 <- foreach(i = 1:B, .packages = c("tidyverse", "secr", "survival"), .combine = c) %dopar% {
  pop <- prune(population(30, 36500), t0)
  
  survhist <- survey_sim(x0, y0, x1, y1, pop, times, t0)$survhist
  coef(survreg(survhist ~ 1, dist = "exponential"))
}


stopCluster(cl)
Sys.time()