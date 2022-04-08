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


cl <- makeCluster(20)
registerDoParallel(cl)
getDoParWorkers()

# Simulate populations of different abundances - with surveys
set.seed(012345)
Sys.time()

N_expected <- 1:200
B <- length(N_expected)

output <- foreach(i = 1:B, .packages = c("tidyverse", "secr")) %dopar% {
  N_exp <- N_expected[i]
  pop <- prune(population(N_exp, 365 * 101), 36500 - 30)
  
  capthist <- survey_sim(x0, y0, x1, y1, pop, times = times, HHN, pars, t0=36500-30)
  mod <- secr::secr.fit(capthist = capthist$secr_capthist, mask = mask)
  
  list(expected = N_exp,
       present = pres_det(pop, 36500)$pres,
       detectable = pres_det(pop, 36500)$det,
       estimated = region.N(mod)[2, 1]
  )
}

df <- data.frame(expected = sapply(output, function(item) item$expected),
                 present = sapply(output, function(item) item$present),
                 detectable = sapply(output, function(item) item$detectable),
                 estimated = sapply(output, function(item) item$estimated))

saveRDS(df, "exp_pres_det_est.rds")
Sys.time()
beep()


# Simulate populations with different mean dropping life time
print("Simulation 2...")
mean_lifes <- 1:365
B <- length(mean_lifes)

output <- foreach(i = 1:B,  .packages = c("tidyverse", "secr")) %dopar% {
  mean_life <- mean_lifes[i]
  pop <- population(30, 365 * 101, mean_life = mean_life)
  
  detectable_times <- sapply(pop,
                             function(ind) max(ind$degrade_times) - min(ind$drop_times))
  
  mean_detectable_time <- mean(detectable_times)
  se <- sd(detectable_times) / sqrt(length(detectable_times))
  
  slope <- mean_detectable_time / (MEAN_STAY + (1 / DROP_RATE))
  
  capthist <- survey_sim(x0, y0, x1, y1, pop, times = times, HHN, pars, t0=36500-30)
  mod <- secr::secr.fit(capthist = capthist$secr_capthist, mask = mask)
  
  list(estimated = region.N(mod)[2, 1],
       present = pres_det(pop, 36500)$pres,
       detectable = pres_det(pop, 36500)$pres,
       mean_life = mean_life,
       mean_detectable_time = mean_detectable_time,
       se = se,
       slope = slope
  )
}

df <- data.frame(mean_life = sapply(output, function(item) item$mean_life),
                 mean_detectable_time = sapply(output, function(item) item$mean_detectable_time),
                 se = sapply(output, function(item) item$se),
                 slope = sapply(output, function(item) item$slope),
                 present = sapply(output, function(item) item$present),
                 estimated = sapply(output, function(item) item$estimated),
                 detectable = sapply(output, function(item) item$detectable))

saveRDS(df, "dropping_life.rds")
Sys.time()



# Simulate populations with different dropping rates
print("Simulation 3...")
dropping_rate <- seq(1/7, 2, length.out = 200)
B <- length(dropping_rate)

output <- foreach(i = 1:B,  .packages = c("tidyverse", "secr")) %dopar% {
  drop_rate <- dropping_rate[i]
  pop <- population(30, 365 * 101, drop_rate = drop_rate)
  
  detectable_times <- sapply(pop,
                             function(ind) max(ind$degrade_times) - min(ind$drop_times))
  
  mean_detectable_time <- mean(detectable_times)
  se <- sd(detectable_times) / sqrt(length(detectable_times))
  
  slope <- mean_detectable_time / (MEAN_STAY + (1 / drop_rate))
  
  capthist <- survey_sim(x0, y0, x1, y1, pop, times = times, HHN, pars, t0=36500-30)
  mod <- secr::secr.fit(capthist = capthist$secr_capthist, mask = mask)
  
  list(drop_rate = drop_rate,
       mean_detectable_time = mean_detectable_time,
       se = se,
       slope = slope,
       estimated = region.N(mod)[2, 1],
       present = pres_det(pop, 36500)$pres,
       detectable = pres_det(pop, 36500)$det
  )
}

df <- data.frame(drop_rate = sapply(output, function(item) item$drop_rate),
                 mean_detectable_time = sapply(output, function(item) item$mean_detectable_time),
                 se = sapply(output, function(item) item$se),
                 slope = sapply(output, function(item) item$slope),
                 estimated = sapply(output, function(item) item$estimated),
                 present = sapply(output, function(item) item$present),
                 detectable = sapply(output, function(item) item$detectable))

saveRDS(df, "drop_rate.rds")
Sys.time()



# Simulate populations with different turnover rates
print("Simulation 4...")
mean_stays <- seq(365, 365 * 2, length.out = 200)
B <- length(mean_stays)

output <- foreach(i = 1:B,  .packages = c("tidyverse", "secr")) %dopar% {
  mean_stay <- mean_stays[i]
  pop <- population(30, 365 * 101, mean_stay = mean_stay)
  
  detectable_times <- sapply(pop,
                             function(ind) max(ind$degrade_times) - min(ind$drop_times))
  
  mean_detectable_time <- mean(detectable_times)
  se <- sd(detectable_times) / sqrt(length(detectable_times))
  
  slope <- mean_detectable_time / (mean_stay + (1 / DROP_RATE))
  
  capthist <- survey_sim(x0, y0, x1, y1, pop, times = times, HHN, pars, t0=36500-30)
  mod <- secr::secr.fit(capthist = capthist$secr_capthist, mask = mask)
  
  list(mean_stay = mean_stay,
       mean_detectable_time = mean_detectable_time,
       se = se,
       slope = slope,
       estimated = region.N(mod)[2, 1],
       present = pres_det(pop, 36500)$pres,
       detectable = pres_det(pop, 36500)$det
  )
}

df <- data.frame(mean_stay = sapply(output, function(item) item$mean_stay),
                 mean_detectable_time = sapply(output, function(item) item$mean_detectable_time),
                 se = sapply(output, function(item) item$se),
                 slope = sapply(output, function(item) item$slope),
                 estimated = sapply(output, function(item) item$estimated),
                 present = sapply(output, function(item) item$present),
                 detectable = sapply(output, function(item) item$detectable))

saveRDS(df, "mean_stay.rds")
stopCluster(cl)
Sys.time()