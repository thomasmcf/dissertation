library(secr)
library(tidyverse)
library(survival)
library(jagsUI)
library(foreach)
library(parallel)
library(doParallel)

# Length of side of study site
SITE_LENGTH <- 200 
# Range parameter for single leopard
RANGE <- 50

# This function simulates the dropping information of one individual leopard, given period of occupation and number of droppings
# Inputs: arrived   - Time that the individual arrived in the study site
#         left      - Time that the individual left the study site
#         droppings - Number of droppings created by the individual during their time in the study area
#                id - Individual identifier
#
# Outputs: A list of six variables:
#          activity_centre - The geographic centre of the individuals activity
#          present         - The time interval the individual was present for
#          drop_times      - The times a dropping was created by the individual
#          surv_times      - The "survival" time of the droppings, which is the time the dropping lasts until degrading beyond genotyping
#          decay_times   - The time at which the dropping degrades beyond genotyping
#                       id - Individual identifier
individual <- function(arrived, left, mean_life, drop_rate, id){
  # Simulate the dropping, survival, and decay times as defined above
  droppings <- rpois(1, drop_rate * (left - arrived))
  drop_times <- runif(droppings, arrived, left)
  surv_times <- rexp(droppings, 1 / mean_life)
  decay_times <- drop_times + surv_times

  # Simulate the activity centre of the individual
  # The coordinates are modelled as independent uniform random variables, centred on (0, 0)
  activity_centre <- runif(2, -SITE_LENGTH/2, SITE_LENGTH/2)
  # Simulate the location of the droppings
  # The coordinates are independent normal random variables centered on the activity centre
  location_x <- rnorm(droppings, activity_centre[1], sqrt(RANGE))
  location_y <- rnorm(droppings, activity_centre[2], sqrt(RANGE))
  
  location <- rbind(location_x, location_y)
  # Generate an ID for each dropping, used to uniquely identify them
  drop_id <- paste(id, 1:droppings, sep = "-")
  
  # Return the simulated information in a list
  return(list(activity_centre = activity_centre,
              present = c(arrived, left),
              drop_times = drop_times,
              surv_times = surv_times,
              decay_times = decay_times,
              location = location,
              id = id,
              drop_id = drop_id))
}


# This function simulates the population with a queue (G/G/infinity I think)
# Inputs:  N - Asymptotic expected number of individuals in the population
#          T - Time to simulate the population for
# Outputs: A list of four elements
#          individuals - A data frame with the time interval each individual was present for
#                times - A vector with times at which the population changes
#           population - Each element is the population in the area as of the corresponding time in times
queue <- function(N, mean_stay, T, var_stay){
  # Set parameters of the gamma distribution by matching moments
  shape <- mean_stay ** 2 / var_stay
  rate <-  mean_stay / var_stay
  
  # Initialise the individuals data frame
  individuals <- data.frame(arrival = double(),
                            departure = double())
  
  # Generate the time at which the next individual moves into the area
  # The shape parameter is adjusted in order to specify the stationary expectation
  latest <- rgamma(1, shape / N, rate)
  # Simulate arrival times until the simulation period has elapsed
  while(latest < T){
    # Draw the new arrival's length of stay, and append this to the data frame
    stay_time <- rgamma(1, shape, rate)
    individuals <- rbind(individuals,
                         list(arrival = latest, departure = latest + stay_time))
    
    # Draw next arrival time
    latest <- latest + rgamma(1, shape / N, rate)
  }
  
  # Find how many individuals have been present in total
  n_individs <- dim(individuals)[1]
  # Convert the data frame to a matrix
  individuals <- as.matrix(individuals)
  
  # "Unroll" all times into one vector
  # The first n_individs are the arrival times, the rest are departure times
  times <- as.vector(individuals)
  # Specify the instantaneous change in population at the times in times
  # For the arrival times this is +1, and the departures it is -1
  change <- rep(c(1, -1), c(n_individs, n_individs))
  
  # Order both change and times in chronological order
  change <- change[order(times)]
  times <- sort(times)
  
  # Calculate the abundance at each time point
  population <- cumsum(change)
  
  # Return the simulated information
  list(individuals = individuals,
       times = times,
       change = change,
       population = population)
}


# This function simulates the dropping information of a population of leopards
# Inputs: N - The asymptotic expected number of individuals
#         T - The time at which the SCR study commences
# Outputs: individuals - A list of N individuals with the information generated by individual()
population <- function(N, T, mean_life = 365/2, drop_rate = 0.5, mean_stay = 365, var_stay = (56 / 1.96) ** 2){
  # Simulate occupancy history of the area
  history <- queue(N, mean_stay, T, var_stay)$individuals
  N_tot <- dim(history)[1]
  
  # Initialise a list to store individuals
  individuals <- vector("list", N_tot)
  
  for(i in 1:N_tot){
    present <- history[i, ]
    # Simulate information on each individual
    individuals[[i]] <- individual(present[1], present[2], mean_life, drop_rate, id = i)
  }
  
  return(individuals)
}

# Determine number of individuals present and detectable in population at time
pres_det <- function(population, time){
  N <- length(population)
  N_present <- 0
  N_detectable <- 0
  
  for(i in 1:N){
    ind <- population[[i]]
    # If animal j was present at tyears, increment the value in the presence vector
    if(ind$present[1] <= time & time <= ind$present[2]){
      N_present <- N_present + 1
    }
    
    # If any of animals j's droppings were present at t=100 years, increment the value in the detectable vector
    if(any(ind$drop_times <= time & time <= ind$decay_times)){
      N_detectable <- N_detectable + 1
    }
  }
  
  list("present" = N_present, "detectable" = N_detectable)
}

prune <- function(population, time){
  N <- length(population)
  
  keep <- NULL
  for(i in 1:N){
    ind <- population[[i]]
    
    if(ind$present[1] <= time & time <= ind$present[2]){
      keep <- c(keep, i)
    }
    else if(any(ind$drop_times <= time & ind$decay_times >= time)){
      keep <- c(keep, i)
    }
  }
  
  new_pop <- vector("list", length(keep))
  for(i in 1:length(keep)){
    new_pop[[i]] <- population[[keep[i]]]
  }
  
  return(new_pop)
}

# This function calculates the minimum distance from the point (a, b) to the line segment joining (x0, y0) and (x1, y1)
min_segment_dist <- function(a, b, x0, y0, x1, y1){
  # Calculate the parameter for the point on the segment closest to (a, b)
  t <- - ((x0 - a)*(x1 - x0) + (y0-b)*(y1-y0)) / ((x1-x0)**2 + (y1-y0)**2)
  
  # If t <= 0, (a, b) is closest to (x0, y0)
  if(t <= 0){
    x_min <- x0
    y_min <- y0
  }
  # If t >= 1, it is closest to (x1, y1)
  else if(t >= 1){
    x_min <- x1
    y_min <- y1
  }
  # If 0 < t < 1, it is closest to a point between (x0, y0) and (x1, y1)
  else{
    x_min <- x0 + (x1 - x0) * t
    y_min <- y0 + (y1 - y0) * t
  }
  
  return(list(dist = sqrt((x_min - a) ** 2 + (y_min - b) ** 2),
              x_min = x_min,
              y_min = y_min))
}

# Simulate the capture history for one occasion
#          x - X coordinates of spline points
#          y - Y coordinates of spline points
# population - Simulated Population object
#       time - The tim at which visit takes place
#   occasion - Which occasion is being simulated
#    session - Which session is being simulated
occasion_sim <- function(x0, y0, x1, y1, population, time, occasion, session = 1, t0){
  # The number of individuals in the population object
  N <- length(population)
  transects <- length(x0)
  
  # Empty data frame to hold the capture history for secr
  capthist <- data.frame(Session = integer(),
                         AnimalID = integer(),
                         Occasion = integer(),
                         X = double(),
                         Y = double())
  
  # Empty data frame to hold capture history for survival analysis
  survhist <- data.frame(droppingID = character(),
                         time = integer(),
                         type = character(),
                         present = character())
  
  # Iterate over the population
  for(i in 1:N){
    # Extract individual i
    ind <- population[[i]]
    # Find the number of droppings produced by this individual
    droppings <- length(ind$drop_times)
    
    for(j in 1:transects){
      # Iterate over each dropping
      for(k in 1:droppings){
        # Check if the dropping was detectable at the time of the survey
        if(ind$drop_times[k] <= time){
          # Find the length between the dropping and the transect
          dist <- min_segment_dist(ind$location[1, k], ind$location[2, k], x0[j], y0[j], x1[j], y1[j])
          # Calculate the probability of detection
          prob <- HHN(dist$dist)
          # Simulate the detection/non-detection (Bernoulli)
          detect <- rbinom(1, 1, prob)
          
          # If the dropping was detected record this
          if(detect == 1){
            # If the dropping had not yet degraded, it is recorded for secr
            if(time <= ind$decay_times[k]){
              capthist <- rbind(capthist,
                                list(Session = session,
                                     AnimalID = i,
                                     Occasion = occasion,
                                     X = dist$x_min,
                                     Y = dist$y_min))
            }
            
            survhist <- rbind(survhist,
                              list(droppingID = ind$drop_id[k],
                                   time = time,
                                   type = ifelse(time <= ind$decay_times[k], "detected", "degraded"),
                                   present = ifelse(ind$drop_times[k] <= t0, TRUE, FALSE)))
          }
        }
      }
    }
  }
  
  survhist <- distinct(survhist)
  
  return(list(capthist = capthist,
              survhist = survhist))
}


# Simulate survey
survey_sim <- function(x0, y0, x1, y1, population, times, t0){
  n_occasions <- length(times)
  occasion <- occasion_sim(x0, y0, x1, y1, population, times[1], occasion = 1, t0 = t0)
  
  capthist <- occasion$capthist
  survhist <- occasion$survhist
  
  for(t in 2:n_occasions){
    occasion <- occasion_sim(x0, y0, x1, y1, population, times[t], occasion = t, t0 = t0)
    
    capthist <- rbind(capthist,
                      occasion$capthist)
    
    survhist <- rbind(survhist,
                      occasion$survhist)
    
  }
  
  captured_ids <- unique(capthist$AnimalID)
  ids <- cbind(sort(captured_ids), 1:length(captured_ids))
  
  for(i in 1:nrow(capthist)){
    old_id <- capthist[i, "AnimalID"]
    new_id <- ids[which(ids[, 1] == old_id), 2]
    
    capthist[i, "AnimalID"] <- new_id
  }
  
  capthist <- make.capthist(capthist, transect_list, fmt="XY")
  survhist <- make.surv(survhist, t0, times)
  
  return(list(capthist = capthist,
              survhist = survhist))
}


make.surv <- function(survhist, t0, times){
  survhist <- survhist %>%
    pivot_wider(names_from = time, values_from = type, values_fill = "no detection")
  
  survhist <- survhist %>%
    mutate(degradation_observed = apply(survhist, 1, function(row) "degraded" %in% row),
           detections = apply(survhist, 1, function(row) sum("detected" == row))) %>%
    filter(!(present & detections == 0))%>%
    filter(!(!present & !degradation_observed & (detections == 1)))
  
  droppings <- nrow(survhist)
  
  time <- numeric(droppings)
  time2 <- numeric(droppings)
  event <- numeric(droppings)
  
  for(i in 1:droppings){
    if(survhist$present[i]){
      t_last <- times[max(which(survhist[i, ] == "detected")) - 2]
      
      time[i] <- t_last - t0
      time2[i] <- NA
      
      event[i] <- 0
    }
    
    else if(survhist$degradation_observed[i]){
      t_deg <- times[min(which(survhist[i, ] == "degraded")) - 2]
      
      if(survhist$detections[i] == 0){
        time[i] <- t_deg - t0
        time2[i] <- NA
        
        event[i] <- 2
      }
      
      else{
        t_init <- times[min(which(survhist[i, ] == "detected")) - 2]
        t_last <- times[max(which(survhist[i, ] == "detected")) - 2]
        
        if(t_init == t_last){
          time[i] <- t_deg
          time2[i] <- NA
          
          event[i] <- 2
        }
        
        else{
          time[i] <- t_last - t_init
          time2[i] <- t_deg - t0
          
          event[i] <- 3 
        }
      }
    }
    
    else{
      t_init <- times[min(which(survhist[i, ] == "detected")) - 2]
      t_last <- times[max(which(survhist[i, ] == "detected")) - 2]
      
      time[i] <- t_last - t_init
      time2[i] <- NA
      
      event[i] <- 0
    }
  }
  
  Surv(time = time, time2 = time2, event = event, type = "interval")
}

HPDI <- function(x, p){
  N <- length(x)
  n <- round(N * p)
  
  x <- sort(x)
  
  min_i <- 1
  min_width <- x[n + 1] - x[1]
  
  for(i in 2:(N-n)){
    lower <- x[i]
    upper <- x[i + n]
    
    width <- upper - lower
    
    if(width < min_width){
      min_width <- width
      min_i <- i
    }
  }
  
  return(x[c(min_i, min_i + n)])
}


lambda <- function(d, lambda0, sigma){
  lambda0 * exp(-d**2 / (2 * sigma**2))
}

HHN <- function(d, lambda0 = 1, sigma = 1){
  1 - exp(-lambda(d, lambda0, sigma))
}