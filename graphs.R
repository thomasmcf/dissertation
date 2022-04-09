library(tidyverse)
library(gridExtra)

source("simulation functions.R")

set.seed(12345)
pop <- population(30, 365 * 101, mean_stay = 365, drop_rate = 0.5, mean_life = 365 / 2)
N <- length(pop)

# Create a blank plot - this will eventually show the spatial locations of activity centres and droppings
plot(c(-SITE_LENGTH/2, SITE_LENGTH/2, -SITE_LENGTH/2, SITE_LENGTH/2),
     c(-SITE_LENGTH/2, SITE_LENGTH/2, SITE_LENGTH/2, -SITE_LENGTH/2),
     type = 'n', xlab = "Longitude", ylab = "Latitude")

for(i in 1:N){
  ind <- pop[[i]]
  # Extract the droppings which are still present at the commencation of the SCR study
  locations <- matrix(ind$location[, which(ind$drop_times <= 100 * 365 & ind$degrade_times >= 100 * 365)], nrow = 2)
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


# EXP PRES DET EST
df <- readRDS("exp_pres_det_est.rds")

detectable_times <- sapply(pop,
                           function(ind) max(ind$degrade_times) - min(ind$drop_times))

mean_detectable_time <- mean(detectable_times)
slope <- mean_detectable_time / (365 + 2)

plot_df <- with(df, data.frame(Y = c(detectable, estimated),
                               Type = rep(c("Detectable", "Estimated"), c(200, 200)),
                               X = rep(present, 2)))

ggplot(plot_df) + geom_point(aes(x = X, y = Y, col = Type), size = 2.25) +
  geom_segment(aes(x = 0, xend = max(X), y = 0, yend = max(X)), size = 1) +
  geom_segment(aes(x = 0, xend = max(X), y = 0, yend = slope * max(X)), size = 1, linetype = 2) +
  xlab("Leopards Present") + 
  ylab("Leopards Estimated") + 
  theme_bw() +
  scale_x_continuous(expand = c(0.025, 0.025)) +
  scale_y_continuous(expand = c(0.025, 0.025)) +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        panel.border = element_blank(),
        axis.line = element_line(colour = 'black'))

# Mean dropping life time
df <- readRDS("dropping_life.rds")

p1 <- ggplot(df) + geom_point(aes(x = mean_life, y = mean_detectable_time)) +
  xlab("Mean Dropping Lifetime") +
  ylab("Mean Time Detectable") +
  theme_bw() +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        panel.border = element_blank(),
        axis.line = element_line(colour = 'black'))

p2 <- ggplot(df) + geom_point(aes(x = mean_life, y = estimated)) +
  xlab("Mean Dropping Lifetime") +
  ylab("Estimated Abundance") +
  theme_bw() +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        panel.border = element_blank(),
        axis.line = element_line(colour = 'black'))

grid.arrange(p1, p2)

# Dropping deposition rate
df <- readRDS("drop_rate.rds")

p1 <- ggplot(df) + geom_point(aes(x = drop_rate, y = mean_detectable_time)) +
  xlab("Dropping Deposition Rate") +
  ylab("Mean Time Detectable") +
  theme_bw() +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        panel.border = element_blank(),
        axis.line = element_line(colour = 'black'))

p2 <- ggplot(df) + geom_point(aes(x = drop_rate, y = estimated)) +
  xlab("Dropping Deposition Rate") +
  ylab("Estimated Abundance") +
  theme_bw() +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        panel.border = element_blank(),
        axis.line = element_line(colour = 'black'))

grid.arrange(p1, p2)

# Mean Stay Time
df <- readRDS("mean_stay.rds")

p1 <- ggplot(df) + geom_point(aes(x = mean_stay, y = mean_detectable_time)) +
  xlab("Mean Stay") +
  ylab("Mean Time Detectable") +
  theme_bw() +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        panel.border = element_blank(),
        axis.line = element_line(colour = 'black'))

p2 <- ggplot(df) + geom_point(aes(x = mean_stay, y = estimated)) +
  xlab("Mean Stay") +
  ylab("Estimated Abundance") +
  theme_bw() +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        panel.border = element_blank(),
        axis.line = element_line(colour = 'black'))

grid.arrange(p1, p2)