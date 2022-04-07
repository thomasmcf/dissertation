library(ggplot2)

set.seed(12345)
pop <- population(30, 365 * 101)
N <- length(pop)

# Create a blank plot - this will eventually show the spatial locations of activity centres and droppings
plot(c(-SITE_LENGTH/2, SITE_LENGTH/2, -SITE_LENGTH/2, SITE_LENGTH/2),
     c(-SITE_LENGTH/2, SITE_LENGTH/2, SITE_LENGTH/2, -SITE_LENGTH/2),
     type = 'n', xlab = "Longitude", ylab = "Latitude")

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

plot <- ggplot() +
  geom_point(aes(x = N_present, y = N_detectable)) +
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