source("simulation functions.R")

sims <- read.csv("simulated.csv")

# Simulate the study site for three years, with N expected individuals
pop <- population(30, 365 * 101)

# Create a blank plot - this will eventually show the spatial locations of activity centres and droppings
plot(c(-SITE_LENGTH/2, SITE_LENGTH/2, -SITE_LENGTH/2, SITE_LENGTH/2),
     c(-SITE_LENGTH/2, SITE_LENGTH/2, SITE_LENGTH/2, -SITE_LENGTH/2),
     type = 'n', xlab = "Longitude", ylab = "Latitude")

# Iterate over the population
set.seed(12345)

for(i in 1:N){
  ind <- pop[[i]]
  # Extract the droppings which are still present at the commencation of the SCR study
  locations <- matrix(ind$location[, which(ind$drop_times <= 100 * 365 & ind$decay_times >= 100 * 365)], nrow = 2)
  # Plot these locations
  points(locations[1,], locations[2,], pch = 19, cex = 0.5, col = 'blue')

  # Plot the activity centres, green if it's present at the start of the study, red if it is not
  if(ind$present[1] <= 100 * 365 & ind$present[2] >= 100 * 365){
    points(ind$activity_centre[1], ind$activity_centre[2], pch = 19, cex = 1.5, col = 'green')
  }
  else if(length(locations) > 0){
    points(ind$activity_centre[1], ind$activity_centre[2], pch = 19, cex = 1.5, col = 'red')
  }
}

# Add a legend
legend(-100, 100, c("Present", "Departed", "Dropping"),
       col = c('green', 'red', 'blue'),
       pch = 19)

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


cl <- makeCluster(2)
registerDoParallel(cl)
getDoParWorkers()


# Simulate populations of different abundances - with surveys
N_expected <- 1:200
B <- length(N_expected)

print("Simulation 1...")
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
                 estimated = sapply(output, function(item) item$esimated))

saveRDS(df, "exp_pres_det_est.rds")



# Simulate populations with different mean dropping life time
print("Simulation 2...")
mean_lifes <- seq(1, 365, 365)
B <- length(mean_life)

output <- foreach(i = 1:B, .packages = "tidyverse") %dopar% {
  mean_life <- mean_lifes[i]
  pop <- population(N_exp, 365 * 101, mean_life = mean_life)
  
  detectable_times <- sapply(pop,
                                  function(ind) max(ind$degrade_times) - min(ind$drop_times))
  
  mean_detectable_time <- mean(detectable_times)
  se <- sd(detectable_times) / sqrt(length(detectable_times))
  
  slope <- mean_detectable_time / (MEAN_STAY + (1 / DROP_RATE))
  
  list(mean_life = mean_life,
       mean_detectable_time = mean_detectable_time,
       se = se,
       slope = slope
  )
}

df <- data.frame(mean_life = sapply(output, function(item) item$mean_life),
                 mean_detectable_time = sapply(output, function(item) item$mean_detectable_time),
                 se = sapply(output, function(item) item$se),
                 slope = sapply(output, function(item) item$slope))

saveRDS(df, "dropping_life.rds")



# Simulate populations with different dropping rates
print("Simulation 3...")
dropping_rate <- seq(1/7, 2, 50)
B <- length(dropping_rate)

output <- foreach(i = 1:B, .packages = "tidyverse") %dopar% {
  drop_rate <- dropping_rate[i]
  pop <- population(N_exp, 365 * 101, drop_rate = drop_rate)
  
  detectable_times <- sapply(pop,
                             function(ind) max(ind$degrade_times) - min(ind$drop_times))
  
  mean_detectable_time <- mean(detectable_times)
  se <- sd(detectable_times) / sqrt(length(detectable_times))
  
  slope <- mean_detectable_time / (MEAN_STAY + (1 / drop_rate))
  
  list(drop_rate = drop_rate,
       mean_detectable_time = mean_detectable_time,
       se = se,
       slope = slope
  )
}

df <- data.frame(drop_rate = sapply(output, function(item) item$drop_rate),
                 mean_detectable_time = sapply(output, function(item) item$mean_detectable_time),
                 se = sapply(output, function(item) item$se),
                 slope = sapply(output, function(item) item$slope))

saveRDS(df, "drop_rate.rds")



# Simulate populations with different turnover rates
print("Simulation 4...")
mean_stays <- seq(30, 365 * 2, 50)
B <- length(mean_stays)

output <- foreach(i = 1:B, .packages = "tidyverse") %dopar% {
  mean_stay <- mean_stays[i]
  pop <- population(N_exp, 365 * 101, mean_stay = mean_stay)
  
  detectable_times <- sapply(pop,
                             function(ind) max(ind$degrade_times) - min(ind$drop_times))
  
  mean_detectable_time <- mean(detectable_times)
  se <- sd(detectable_times) / sqrt(length(detectable_times))
  
  slope <- mean_detectable_time / (mean_stay + (1 / DROP_RATE))
  
  list(mean_stay = mean_stay,
       mean_detectable_time = mean_detectable_time,
       se = se,
       slope = slope
  )
}

df <- data.frame(mean_stay = sapply(output, function(item) item$mean_stay),
                 mean_detectable_time = sapply(output, function(item) item$mean_detectable_time),
                 se = sapply(output, function(item) item$se),
                 slope = sapply(output, function(item) item$slope))

saveRDS(df, "mean_stay.rds")