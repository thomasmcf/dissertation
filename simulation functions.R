library(secr)
library(MASS)
library(tidyverse)
library(beepr)

setwd("C:/Users/thmcf/Documents/Programming/Maths/Dissertation/dissertation/")

# Mean life time of scat
MEAN_LIFE <- 365 / 2
# Mean droppings left per day
DROP_RATE <- 0.5
# Mean time before leaving area
MEAN_STAY <- 365
# Variance in the time spent in area
VAR_STAY <- (56 / 1.96) ** 2
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
#          degrade_times   - The time at which the dropping degrades beyond genotyping
#          vis_times       - The times which droppings are visible for after degrading beyond genotyping
#          decay_times     - The times at which droppings decay and are no longer detectable
#                       id - Individual identifier
individual <- function(arrived, left, id){
  # Simulate the dropping, survival, and decay times as defined above
  droppings <- rpois(1, DROP_RATE * (left - arrived))
  drop_times <- runif(droppings, arrived, left)
  surv_times <- rexp(droppings, 1 / MEAN_LIFE)
  degrade_times <- drop_times + surv_times
  vis_times <- rexp(droppings, 1 / MEAN_LIFE)
  decay_times <- degrade_times + vis_times
  
  # Simulate the activity centre of the individual
  # The coordinates are modelled as independent uniform random variables, centred on (0, 0)
  activity_centre <- runif(2, -SITE_LENGTH/2, SITE_LENGTH/2)
  # Simulate the location of the droppings
  # The coordinates are independent normal random variables centered on the activity centre
  location <- t(mvrnorm(droppings, activity_centre, RANGE * diag(c(1, 1))))
  # Generate an ID for each dropping, used to uniquely identify them
  drop_id <- paste(id, 1:droppings, sep = "-")
  
  # Return the simulated information in a list
  return(list(activity_centre = activity_centre,
              present = c(arrived, left),
              drop_times = drop_times,
              surv_times = surv_times,
              degrade_times = degrade_times,
              vis_times = vis_times,
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
queue <- function(N, T){
  # Set parameters of the gamma distribution by matching moments
  shape <- MEAN_STAY ** 2 / VAR_STAY
  rate <-  MEAN_STAY / VAR_STAY
  
  # Initialise the individuals data frame
  individuals <- data.frame(arrival = double(),
                            departure = double())
  
  # Generate the time at which the next individual moves into the area
  # The shape paramater is adjusted in order to specify the stationary expectation
  latest <- rgamma(1, shape / N, rate)
  # Simulate arrival times until the simlation period has elapsed
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
population <- function(N, T){
  # Simulate occupancy history of the area
  history <- queue(N, T)$individuals
  N_tot <- dim(history)[1]
  
  # Initialise a list to store individuals
  individuals <- vector("list", N_tot)
  
  for(i in 1:N_tot){
    present <- history[i, ]
    # Simulate information on each individual
    individuals[[i]] <- individual(present[1], present[2], id = i)
  }
  
  return(individuals)
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

# This function calculates the minimum distance from the point (a, b) to a curve taking the form of
# a sequence of linear segments
# x, y are vectors of the points where linear segments start and end
# E.g. segment one runs from (x[1], y[1]) to (x[2], y[2]), segment two from (x[2], y[2]) to (x[3], y[3]), etc.
min_curve_dist <- function(a, b, x, y){
  points <- length(x)
  
  # Calculate distance from point to first segment
  min_dist <- min_segment_dist(a, b, x[1], y[1], x[2], y[2])
  
  for(i in 2:(points - 1)){
    # Calculate distance from point to segment i
    dist_i <- min_segment_dist(a, b, x[i], y[i], x[i+1], y[i+1])
    
    # If this distance is less than the current running minimum, update the minimum
    if(dist_i$dist < min_dist$dist){
      min_dist <- dist_i
    }
  }
  
  return(min_dist)
}

# Simulate the capture history for one occasion
#          x - X coordinates of spline points
#          y - Y coordinates of spline points
# population - Simulated Population object
#       time - The tim at which visit takes place
#   detectfn - Detection function to use
#  detectpar - Parameters of the detection function
#   occasion - Which occasion is being simulated
#    session - Which session is being simulated
occasion_sim <- function(x, y, population, time, detectfn, detectpar, occasion, session = 1){
  # The number of individuals in the population object
  N <- length(population)
  # Empty data frame to hold the capture history for secr
  history <- data.frame(Session = integer(),
                        AnimalID = integer(),
                        Occasion = integer(),
                        X = double(),
                        Y = double())
  
  # Empty data frame to hold capture history for survival analysis
  detections <- data.frame(droppingID = character(),
                           time = integer(),
                           type = character())
  
  # Iterate over the population
  for(i in 1:N){
    # Extract individual i
    ind <- population[[i]]
    # Find the number of droppings produced by this individual
    droppings <- length(ind$drop_times)
    
    # Iterate over each dropping
    for(j in 1:droppings){
      # Check if the dropping was detectable at the time of the survey
      if(ind$drop_times[j] <= time & time <= ind$decay_times[j]){
        # Find the length between the dropping and the transect
        dist <- min_curve_dist(ind$location[1, j], ind$location[2, j], x, y)
        # Calculate the probability of detection
        prob <- detectfn(dist$dist, detectpar)
        # Simulate the detection/non-detection (Bernoulli)
        detect <- rbinom(1, 1, prob)
        
        # If the dropping was detected record this
        if(detect == 1){
          # If the dropping had not yet degraded, it is recorded for secr and survival analysis
          if(time <= ind$degrade_times[j]){
            history <- rbind(history,
                             list(Session = session,
                                  AnimalID = i,
                                  Occasion = occasion,
                                  X = dist$x_min,
                                  Y = dist$y_min))
            
            # And the detection history as an identifiable sample
            detections <- rbind(detections,
                                list(droppingID = ind$drop_id[j],
                                     time = time,
                                     type = "detected"))
          }
          
          # If the dropping was degraded but not decayed, it is only recorded for survival
          else{
            detections <- rbind(detections,
                                list(droppingID = ind$drop_id[j],
                                     time = time,
                                     type = "degraded"))
          }
        }
      }
    }
  }
  
  return(list(capthist = history,
              detections = detections))
}

# Simulate survey
survey_sim <- function(x, y, population, times, detectfn, detectpar){
  n_occasions <- length(times)
  occ <- occasion_sim(x, y, pop, times[1], HHN, pars, occasion = 1)
  
  capthist <- occ$capthist
  detections <- occ$detections
  
  for(t in 2:n_occasions){
    occ <- occasion_sim(x, y, pop, times[t], HHN, pars, t)
    
    capthist <- rbind(capthist,
                      occ$capthist)
    
    detections <- rbind(detections,
                        occ$detections)
  }
  
  ids <- cbind(sort(unique(capthist$AnimalID)), 1:length(unique(capthist$AnimalID)))
  
  for(i in 1:nrow(capthist)){
    capthist[i, "AnimalID"] <- ids[which(ids[, 1] == capthist[i, "AnimalID"]), 2]
  }
  
  secr_capthist <- make.capthist(capthist, transect_list, fmt="XY")
  
  detection <- detections %>%
    pivot_wider(id_cols = time, names_from = droppingID, values_from = type, values_fill = "not-detected") %>%
    t
  
  return(list(secr_capthist = secr_capthist,
              dethist = detection))
}