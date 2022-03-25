setwd("C:/Users/thmcf/Documents/Programming/Maths/Dissertation/dissertation")
source("simulation functions.R")

# Simulate the study site for three years, with N expected individuals
N_exp <- 30
pop <- population(N_exp, 365 * 101)

# Create a blank plot - this will eventually show the spatial locations of activity centres and droppings
plot(c(-SITE_LENGTH/2, SITE_LENGTH/2, -SITE_LENGTH/2, SITE_LENGTH/2),
     c(-SITE_LENGTH/2, SITE_LENGTH/2, SITE_LENGTH/2, -SITE_LENGTH/2),
     type = 'n', xlab = "Longitude", ylab = "Latitude")

# Iterate over the population
set.seed(12345)
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
N_estimated <- numeric(B)
set.seed(196)

transect_list <- read.traps(file = "transects.txt", detector = "transect")
x <- transect_list$x
y <- transect_list$y

x0 <- x[seq(1, 22, 2)]
y0 <- y[seq(1, 22, 2)]
x1 <- x[seq(2, 22, 2)]
y1 <- y[seq(2, 22, 2)]

pars <- list(lambda0 = 1, sigma = 30)

mask <- make.mask(traps = transect_list, buffer = 0, spacing = 10)

# Sys.time() # Time how long the simulations take
# Run simulations - Takes ~7 minutes on my machine
for(i in 1:B){
  # Pick a number of expected individuals
  N_exp <- 30
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
  
  print(i)
}

capthist <- survey_sim(x0, y0, x1, y1, pop, times = c(36500, 36530, 36560), HHN, pars)
mod <- secr.fit(capthist = capthist, buffer = buffer)

# Plot the results
#plot(N_present, N_detectable)

plot <- ggplot() + geom_point(aes(x = N_present, y = N_detectable)) +
  geom_abline(slope = 1, intercept = 0) +
  xlab("Leopards Present") +
  ylab("Leopards Detectable") +
  theme_bw() +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = 'black'))

abline(0, 1)
beep()



N_present <- numeric(20)
N_estimated <- numeric(20)
# Sys.time() # End timing
for(N_exp in seq(10, 200, 10)){
  print(paste("Simulation", N_exp / 10))
  print("    Simulating population")
  pop <- population(N_exp, 365 * 101)
  
  print("    Simulating survey")
  capthist <- survey_sim(x, y, pop, c(36500, 36530, 36560), HHN, pars)
  print("    Fitting model")
  mod <- secr.fit(capthist = capthist$secr_capthist, trace=F, mask = mask)
  
  N_estimated[N_exp / 10] <-  region.N(mod)[2, 1]
  
  for(j in 1:length(pop)){
    # If animal j was present at t=100 years, increment the value in the presence vector
    if(between(365 * 100, pop[[j]]$present[1], pop[[j]]$present[2])){
      N_present[N_exp / 10] <- N_present[N_exp /10] + 1
    }
  }
}