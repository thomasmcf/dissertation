source("simulation functions.R")

# Simulate the study site for three years, with N expected individuals
N_exp <- 30
pop <- population(N_exp, 365 * 3)

# Create a blank plot - this will eventually show the spatial locations of activity centres and droppings
plot(c(-SITE_LENGTH/2, SITE_LENGTH/2, -SITE_LENGTH/2, SITE_LENGTH/2),
     c(-SITE_LENGTH/2, SITE_LENGTH/2, SITE_LENGTH/2, -SITE_LENGTH/2),
     type = 'n', xlab = "Longitude", ylab = "Latitude")

# Iterate over the population
N <- length(pop)
for(i in 1:N){
  ind <- pop[[i]]
  # Extract the droppings which are still present at the commencation of the SCR study
  locations <- matrix(ind$location[, which(ind$decay_times >= 3 * 365)], nrow = 2)
  # Plot these locations
  points(locations[1,], locations[2,], pch = 19, cex = 0.5, col = 'blue')
  
  # Plot the activity centres, green if it's present at the start of the study, red if it is not
  if(ind$present[2] >= 1095){
    points(ind$activity_centre[1], ind$activity_centre[2], pch = 19, cex = 1.5, col = 'green')
  }
  else{
    points(ind$activity_centre[1], ind$activity_centre[2], pch = 19, cex = 1.5, col = 'red')
  }
}

# Add a legend
legend(-100, 100, c("Present", "Departed", "Dropping"),
       col = c('green', 'red', 'blue'),
       pch = 19)
 
# Now explore the issue of "prolonged detectability"
# Define a number of simulations to run, and define vectors to store results
B <- 1
N_present <- numeric(B)
N_detectable <- numeric(B)

# Sys.time() # Time how long the simulations take
# Run simulations - Takes ~7 minutes on my machine
for(i in 1:B){
  # Pick a number of expected individuals
  N_exp <- sample(5:200, 1)
  # Generate 101 years of data
  # The survey is carried out at t=100 years
  # The first 100 years are for burn in, and to ensure no droppings would be present at t=100 that were dropped before t=0
  # The last year is to avoid the tail behaviour that occurs as the simulation ends
  pop <- population(N_exp, 365 * 101)
  
  # Iterate over all animals that were present
  for(j in 1:length(pop)){
    # If animal j was present at t=100 years, increment the value in the presence vector
    if(between(365 * 100, pop[[j]]$present[1], pop[[j]]$present[2])){
      N_present[i] <- N_present[i] + 1
    }
    
    # If any of animals j's droppings were present at t=100 years, increment the value in the detectable vector
    if(any(pop[[j]]$drop_times <= 365 * 100 & pop[[j]]$degrade_times >= 365 * 100)){
      N_detectable[i] <- N_detectable[i] + 1
    }
  }
}

# Plot the results
plot(N_present, N_detectable)
abline(0, 1)
beep()

# Sys.time() # End timing

lambda <- function(d, lambda0, sigma){
  lambda0 * exp(-d**2 / (2 * sigma**2))
}

HHN <- function(d, pars){
  1 - exp(-lambda(d, pars$lambda0, pars$sigma))
}

transect_list <- read.traps(file = "transects.txt", detector = "transect")

pop <- population(30, 365 * 101)
x <- transect_list$x
y <- transect_list$y
pars <- list(lambda0 = 0.5, sigma = 50)

capthist <- survey_sim(x, y, pop, seq(36500, 36500 + 12 * 30, 30))

mask <- make.mask(traps = transect_list, buffer = 0, spacing = 10)

# Fit model (or at least try)
mod <- secr.fit(capthist = capthist$secr_capthist, trace=T, mask = mask)
region.N(mod)